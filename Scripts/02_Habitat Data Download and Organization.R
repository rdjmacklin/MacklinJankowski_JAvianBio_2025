################################################################################

# Script Title: Downloading, organizing, and importing habitat data.

# Note: Much of the code in this script is taken from Dr. Matthew Strimas-Mackey's
# fantastic guide on the analysis of eBird Data, which can be found here:
# https://cornelllabofornithology.github.io/ebird-best-practices/.

################################################################################

## Install required packages (if not already installed).

librarian::shelf("sf", "gdistance", "raster", "exactextractr", "viridis",
                 "tidyverse", "rnaturalearth", "remotes", "landscapemetrics",
                 "units", "terra", "svMisc","fdetsch/MODIS")

# Resolve namespace conflicts

select <- dplyr::select
map <- purrr::map
projection <- raster::projection

# Set S2 to FALSE to facilitate spatial analysis

sf::sf_use_s2(FALSE)

## Import/download several necessary datasets

# Create a .gpkg to store Natural Earth Data

gpkg_filepath <- file.path("./Data/Raw/Natural-Earth/ne_data.gpkg")

# Download naturalearth shapefiles for: 

# Land border w/o lakes

ne_land <- ne_download(scale = 50, category = "cultural",
                       type = "admin_0_countries_lakes",
                       returnclass = "sf") %>%
  filter(CONTINENT == "North America") %>%
  st_set_precision(1e6) %>%
  st_union()

# Country lines (filtered to North America with st_union)

ne_country_lines <- ne_download(scale = 50, category = "cultural",
                                type = "admin_0_boundary_lines_land",
                                returnclass = "sf") %>% 
  st_geometry()

ne_country_lines <- st_intersects(ne_country_lines, ne_land, sparse = FALSE) %>%
  as.logical() %>%
  {ne_country_lines[.]}

# State and provincial boundaries

ne_state_lines <- ne_download(scale = 50, category = "cultural",
                              type = "admin_1_states_provinces_lines",
                              returnclass = "sf") %>%
  filter(ADM0_A3 %in% c("USA", "CAN")) %>%
  mutate(iso_a2 = recode(ADM0_A3, USA = "US", CAN = "CAN")) %>% 
  select(country = ADM0_NAME, country_code = iso_a2)

# Read in county district lines

ne_counties <- ne_download(scale = 10, category = "cultural",
                           type = "admin_2_counties",
                           returnclass = "sf")
# Write to gpkg

write_sf(ne_land, gpkg_filepath, "ne_land")
write_sf(ne_country_lines, gpkg_filepath, "ne_country_lines")
write_sf(ne_state_lines, gpkg_filepath, "ne_state_lines")
write_sf(ne_counties, gpkg_filepath, "ne_counties")


# Download, filter and import a shapefile for the BCR 5 boundary

if(!(file.exists("./Data/Raw/BCR_Terrestrial/BCR_Terrestrial_master_International.shp"))) {
  
  paste0("https://birdscanada.org/download/gislab/bcr_terrestrial_shape.zip") %>% 
    download.file(destfile = "./Data/Raw/BCR_Terrestrial/bcr.zip")
  
  unzip("./Data/Raw/BCR_Terrestrial/bcr.zip", junkpaths = TRUE, exdir = "./Data/Raw/BCR_Terrestrial")
  
  file.remove("./Data/Raw/BCR_Terrestrial/bcr.zip")
}

bcr <- file.path("./Data/Raw/BCR_Terrestrial/BCR_Terrestrial_master_International.shp") %>% 
  read_sf() %>%
  select(bcr_code = BCR, bcr_name = Label) %>% 
  filter(bcr_code == 5) %>%
  
  # Project to the native modis projection
  st_transform(crs = paste("+proj=sinu +lon_0=0 +x_0=0 +y_0=0",
                           "+a=6371007.181 +b=6371007.181 +units=m +no_defs"))

sf_counties <- ne_counties %>%
  filter(REGION == "CA", NAME %in% c("Napa", "Marin", "Solano", "Alameda", "San Mateo", "San Francisco",
                                     "Contra Costa", "Santa Clara", "Sonoma")) %>%
  st_transform(crs = paste("+proj=sinu +lon_0=0 +x_0=0 +y_0=0",
                           "+a=6371007.181 +b=6371007.181 +units=m +no_defs")) %>%
  st_union()

study_region <- st_union(bcr, sf_counties)

## Create a prediction surface of our habitat variables to predict upon when
## occupancy models are built.

laea_crs <- st_crs("+proj=laea +lat_0=51.5 +lon_0=-128")

# Create a raster template covering the region with 1km and 3km resolution

r_1km <- rast(st_transform(study_region, laea_crs), res = c(1000, 1000))

r_3km <- rast(st_transform(study_region, laea_crs), res = c(3000, 3000))

# Fill the raster with 1s inside the study region

r_1km <- rasterize(st_transform(study_region, laea_crs), r_1km, values = 1) %>% 
  setNames("study_region")

r_3km <- rasterize(st_transform(study_region, laea_crs), r_3km, values = 1) %>% 
  setNames("study_region")

# Save for later use

r_1km <- writeRaster(r_1km, "./Data/Cleaned/Prediction-Grids/prediction-grid_1km.tif",
                     overwrite = TRUE,
                     gdal = "COMPRESS=DEFLATE")

r_3km <- writeRaster(r_3km, "./Data/Cleaned/Prediction-Grids/prediction-grid_3km.tif",
                     overwrite = TRUE,
                     gdal = "COMPRESS=DEFLATE")

## Download MODIS data

# Get list of tiles required to cover this BCR

tiles <- getTile(bcr)

tiles@tile

# Enter LOGIN information for NASA Earthdata account - register at 
# https://urs.earthdata.nasa.gov/users/new

MODIS::EarthdataLogin(usr = "INSERT USERNAME HERE", pwd = "INSERT PASSWORD HERE")

# Define start- and end-date for MODIS data download

begin_year <- "2012.01.01"
end_year <- "2021.12.31"

# Download tiles and combine into a single raster for each year

tifs <- runGdal(product = "MCD12Q1", collection = "061", SDSstring = "01", 
                tileH = tiles@tileH, tileV = tiles@tileV, 
                begin = begin_year, end = end_year, 
                outDirPath = "./Data/Raw", job = "MODIS",
                MODISserverOrder = "LPDAAC") %>%
  pluck("MCD12Q1.061") %>% 
  unlist()

# Reformat names to display year and rename

new_names <- format(as.Date(names(tifs)), "%Y") %>% 
  sprintf("modis_mcd12q1_umd_%s.tif", .) %>% 
  file.path(dirname(tifs), .)

file.rename(tifs, new_names)

# Import into a stacked raster with all years data

landcover <- list.files("./Data/Raw/MODIS", "^modis_mcd12q1_umd", 
                        full.names = TRUE) %>%
  stack()

# Fix strange error in layer naming by referencing the filepath stored in each layer.

for(i in 1:length(names(landcover))) {
  
  names(landcover)[i] <- landcover@layers[[i]]@file@name
  
}

# Label layers with year

landcover <- names(landcover) %>% 
  str_extract("(?<=modis_mcd12q1_umd_)[0-9]{4}") %>% 
  paste0("y", .) %>% 
  setNames(landcover, .)

landcover

# Extract minimum and maximum year with land-cover data available

min_lc_year <- names(landcover) %>% 
  str_extract("[0-9]{4}") %>% 
  as.integer() %>% 
  min()

max_lc_year <- names(landcover) %>% 
  str_extract("[0-9]{4}") %>% 
  as.integer() %>% 
  max()

# Create a vector of landcover classes and names to use later

lc_classes <- as.data.frame(cbind(c(0:10, 12, 13, 15, 255),
                                  c("water", "evergreen_needleleaf",
                                    "evergreen_broadleaf", "deciduous_needleleaf",
                                    "deciduous_broadleaf", "mixed_forest",
                                    "closed_shrubland", "open_shrubland",
                                    "woody_savanna", "savanna", "grassland",
                                    "cropland", "urban", "nonvegetated",
                                    "unclassified")))

names(lc_classes) <- c("class", "label")

lc_classes$class <- as.integer(lc_classes$class)

## Import filtered and zero-filled eBird Data - this can be for either species 
## as both zero-filledfiles contain all checklist locations we want to extract 
## habitat data for.

CBCH_eBird <- read_csv("./Data/Cleaned/eBird-Basic-Dataset/Zero-Filled-and-Filtered/CBCH/ebd_CBCH_BCR5_zf.csv") %>%
  
  # Filter to checklists to be used in analysis: 2012-2021
  filter(year %in% c(min_lc_year:max_lc_year)) %>%
  mutate(year_lc = paste0("y", year))

## Define locations at which to extract landcover data based on locations of 
## eBird checklists from our filtered and zero-filled eBird Basic Dataset

## Create a set of buffered eBird checklist locations

eBird_buff <- CBCH_eBird %>% 
  
  # Select distinct location-year combinations
  distinct(year, locality_id, latitude=round(latitude,4), longitude=round(longitude,4)) %>% 
  
  # Create a variable indicating which landcover year data should be used. This
  # can be modified to allow inclusion of years of eBird data beyond the maximum 
  # year available for landcover data by specifying those years use the most recent
  # year of landcover data available.
  mutate(year_lc = paste0("y", year)) %>% 
  
  # Convert to spatial features
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
  
  # Transform to MODIS projection
  st_transform(crs = projection(landcover)) %>% 
  
  # Buffer to create neighborhood around each point
  st_buffer(dist = set_units(1.5, "km"))

## Extract landcover values around each buffered checklist to be used in
## proportion calculations

lsm <- NULL

for (i in seq_len(nrow(eBird_buff))) {
  
  buffer_i <- eBird_buff[i,]
  
  year <- as.character(buffer_i$year_lc)
  
  # crop and mask landcover raster so all values outside buffer are missing
  
  lsm[[i]] <- crop(landcover[[year]], buffer_i) %>% 
    
    mask(buffer_i) %>% 
    
    # calculate landscape metrics
    calculate_lsm(level = "class", metric = c("pland", "ed")) %>% 
    
    # add variables to uniquely identify each point
    mutate(locality_id = buffer_i$locality_id, 
           year_lc = buffer_i$year_lc)
  
  progress(i, max.value = nrow(eBird_buff))
  
}

lsm <- bind_rows(lsm) %>% 
  select(locality_id, year_lc, class, metric, value)

lsm_wide <- lsm %>% 
  # fill missing classes with zeros
  complete(nesting(locality_id, year_lc),
           class = lc_classes$class,
           metric = c("ed", "pland"),
           fill = list(value = 0)) %>% 
  # bring in more descriptive names
  inner_join(lc_classes, by = "class") %>% 
  # transform from long to wide format
  pivot_wider(values_from = value,
              names_from = c(class, label, metric),
              names_glue = "{metric}_c{str_pad(class, 2, pad = '0')}_{label}",
              names_sort = TRUE) %>% 
  arrange(locality_id, year_lc)

## Create objects to fill with covariate data.

yearly_vars <- lsm_wide

## Now we want to link elevation-related values to each checklist.

# Load rasters

elevation <- rast("./Data/Raw/Elevation/elevation_1KMmd_GMTEDmd.tif")
slope <- rast("./Data/Raw/Elevation/slope_1KMmd_GMTEDmd.tif")
eastness <- rast("./Data/Raw/Elevation/eastness_1KMmd_GMTEDmd.tif")
northness <- rast("./Data/Raw/Elevation/northness_1KMmd_GMTEDmd.tif")

# Transform the buffers to the elevation raster's CRS

eBird_buff <- st_transform(eBird_buff, crs = st_crs(elevation))

# Extract the mean and standard deviation of the elevation in the neighbourhoods.

yearly_vars <- exact_extract(elevation, eBird_buff, fun = c("mean", "stdev"),
                             progress = TRUE) %>% 
  # Add variables to uniquely identify each point
  mutate(locality_id = eBird_buff$locality_id, year_lc = eBird_buff$year_lc) %>%
  select(locality_id, year_lc, elevation = mean, elevation_sd = stdev) %>%
  left_join(yearly_vars, ., by = c("locality_id", "year_lc"))

yearly_vars <- exact_extract(slope, eBird_buff, fun = c("mean"),
                             progress = TRUE) %>%
  as.data.frame() %>%
  setNames("mean") %>% 
  # add variables to uniquely identify each point
  mutate(locality_id = eBird_buff$locality_id, year_lc = eBird_buff$year_lc) %>% 
  select(locality_id, year_lc, slope = mean) %>%
  left_join(yearly_vars, ., by = c("locality_id", "year_lc"))

yearly_vars <- exact_extract(northness, eBird_buff, fun = c("mean"),
                             progress = TRUE) %>%
  as.data.frame() %>%
  setNames("mean") %>% 
  # add variables to uniquely identify each point
  mutate(locality_id = eBird_buff$locality_id, year_lc = eBird_buff$year_lc) %>% 
  select(locality_id, year_lc, northness = mean) %>%
  left_join(yearly_vars, ., by = c("locality_id", "year_lc"))

yearly_vars <- exact_extract(eastness, eBird_buff, fun = c("mean"),
                             progress = TRUE) %>%
  as.data.frame() %>%
  setNames("mean") %>% 
  # add variables to uniquely identify each point
  mutate(locality_id = eBird_buff$locality_id, year_lc = eBird_buff$year_lc) %>% 
  select(locality_id, year_lc, eastness = mean) %>%
  left_join(yearly_vars, ., by = c("locality_id", "year_lc"))

## Now we want to extract the VIIRS radiance value at each checklist location

# Import the VIIRS dataset

VIIRS <- list.files("./Data/Raw/VIIRS", 
                    full.names = TRUE) %>%
  stack()


VIIRS <- names(VIIRS) %>% 
  substr(start = 13, stop = 16) %>% 
  paste0("y", .) %>% 
  setNames(VIIRS, .)

VIIRS

# Create an object to fill with the nighttime light value at each locality

nighttime_light <- eBird_buff[,c("locality_id", "year_lc")] %>%
  st_drop_geometry()

# Create a new column to fill with the nighttime light value

nighttime_light$nighttime_light <- NA

# Transform the buffers to the road density raster's CRS

eBird_buff <- st_transform(eBird_buff, crs = st_crs(VIIRS))

# Cycle through each location, and take the mean nighttime light value in the
# location neighbourhood.

for(i in seq_len(nrow(eBird_buff))) {
  
  nighttime_light$nighttime_light[i] <- exact_extract(VIIRS[[paste0("y", eBird_buff$year[i])]],
                                                      eBird_buff[i,], fun = c("mean"),
                                                      progress = FALSE)
  
  progress(which(seq_len(nrow(eBird_buff)) == i), max.value = nrow(eBird_buff))
  
}

# Add nighttime light to the other environmental variables.

yearly_vars <- left_join(yearly_vars, nighttime_light, by = c("locality_id", "year_lc"))

# Load road vector image.

road_dens <- vect("./Data/Raw/GRIP4-Vector/GRIP4_region1.gdb")

# We will convert this into a raster with length being the value of each cell.
# We want the output units to be m/km^2, so we will do this with our 1km^2
# resolution basemap.

r_1km <- terra::project(r_1km, crs(road_dens))

road_dens <- crop(road_dens, r_1km) 

road_dens <- rasterizeGeom(road_dens, r_1km, fun = "length", unit = "m")

# Transform the buffers to the road density raster's CRS

eBird_buff <- st_transform(eBird_buff, crs = st_crs(road_dens))

# Extract road density form our raster at each neighbourhood.

yearly_vars <- exact_extract(road_dens, eBird_buff, fun = c("mean"),
                             progress = TRUE) %>%
  as.data.frame() %>%
  setNames("mean") %>%
  # add variables to uniquely identify each point
  mutate(locality_id = eBird_buff$locality_id, year_lc = eBird_buff$year_lc) %>% 
  select(locality_id, year_lc, road_dens = mean) %>%
  left_join(yearly_vars, ., by = c("locality_id", "year_lc"))

# Load canopy height raster

canopy_height <- rast("./Data/Raw/GEDI/Forest_height_2019_NAM.tif")

canopy_height <- study_region %>% 
  st_buffer(dist = 10000) %>% 
  st_transform(crs = crs(canopy_height)) %>% 
  crop(canopy_height, .)

canopy_height <- extent(-135, -120.7977, 36.5, 51) %>%
  crop(canopy_height, .) %>%
  classify(., cbind(101,NA))

# Transform the buffers to the canopy height raster's CRS

eBird_buff <- st_transform(eBird_buff, crs = st_crs(canopy_height))

# Extract the canopy height in each neighbourhood.

yearly_vars <- exact_extract(canopy_height, eBird_buff, fun = c("mean"), progress = TRUE) %>%
  as.data.frame() %>%
  setNames("mean") %>%
  # add variables to uniquely identify each point
  mutate(locality_id = eBird_buff$locality_id, year_lc = eBird_buff$year_lc) %>% 
  select(locality_id, year_lc, canopy_height = mean) %>%
  left_join(yearly_vars, ., by = c("locality_id", "year_lc"))


## Now we want to extract the climate values at each checklist location. These 
## will be monthly average so we will recreate our buffers including month
## as a unique factor.

eBird_buff_monthly <- CBCH_eBird %>% 
  
  mutate(month = month(observation_date)) %>%
  
  # Select distinct location-year combinations
  distinct(year, month, locality_id, latitude=round(latitude,4), longitude=round(longitude,4)) %>% 
  
  # Create a variable indicating which landcover year data should be used. This
  # can be modified to allow inclusion of years of eBird data beyond the maximum 
  # year available for landcover data by specifying those years use the most recent
  # year of landcover data available.
  mutate(year_lc = paste0("y", year)) %>% 
  
  # Convert to spatial features
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
  
  # Transform to MODIS projection
  st_transform(crs = projection(landcover)) %>% 
  
  # Buffer to create neighborhood around each point
  st_buffer(dist = set_units(1.5, "km"))

# To make extracting raster values faster, we will break up the buffers into
# groups of buffers that have the same year and month data.

eBird_buff_months <- list()

# List all year-month combinations

i_months <- paste0(rep(sort(unique(eBird_buff_monthly$year)), each = length(unique(eBird_buff_monthly$month))), "_", rep(sort(unique(eBird_buff_monthly$month)), times = length(unique(eBird_buff_monthly$year))))

# Loop through and filter.

for(i in i_months) {
  
  eBird_buff_months[[i]] <- eBird_buff_monthly %>%
    filter(year == as.numeric(substr(i, 1, 4)), month == as.numeric(sub('.*_', '', i)))
  
}

# Import the precipitation dataset

daymet_prcp <- list.files("./Data/Raw/Daymet", 
                          full.names = TRUE,
                          pattern = "daymet_v4_prcp_monttl*") %>%
  stack()


daymet_prcp <- names(daymet_prcp) %>% 
  substr(start = 26, stop = 29) %>% 
  paste0("y", .,".", ifelse(nchar(names(daymet_prcp)) == 33, substr(names(daymet_prcp), 33, 33), substr(names(daymet_prcp), 33, 34))) %>% 
  setNames(daymet_prcp, .)

daymet_prcp

# Import the min temperature dataset

daymet_tmin <- list.files("./Data/Raw/Daymet", 
                          full.names = TRUE,
                          pattern = "daymet_v4_tmin_monavg*") %>%
  stack()


daymet_tmin <- names(daymet_tmin) %>% 
  substr(start = 26, stop = 29) %>% 
  paste0("y", .,".", ifelse(nchar(names(daymet_tmin)) == 33, substr(names(daymet_tmin), 33, 33), substr(names(daymet_tmin), 33, 34))) %>% 
  setNames(daymet_tmin, .)

daymet_tmin


# Import the max temperature dataset

daymet_tmax <- list.files("./Data/Raw/Daymet", 
                          full.names = TRUE,
                          pattern = "daymet_v4_tmax_monavg*") %>%
  stack()


daymet_tmax <- names(daymet_tmax) %>% 
  substr(start = 26, stop = 29) %>% 
  paste0("y", .,".", ifelse(nchar(names(daymet_tmax)) == 33, substr(names(daymet_tmax), 33, 33), substr(names(daymet_tmax), 33, 34))) %>% 
  setNames(daymet_tmax, .)

daymet_tmax

# Loop through and extract from each layer.

for(i in names(eBird_buff_months)) {
  
  eBird_buff_months[[i]] <- st_transform(eBird_buff_months[[i]], crs = st_crs(daymet_prcp))
  
  
  eBird_buff_months[[i]]$prcp <- exact_extract(daymet_prcp[[paste0("y", unique(eBird_buff_months[[i]]$year),".", unique(eBird_buff_months[[i]]$month))]],
                                               eBird_buff_months[[i]], fun = c("mean"),
                                               progress = FALSE)
  
  eBird_buff_months[[i]]$tmin <- exact_extract(daymet_tmin[[paste0("y", unique(eBird_buff_months[[i]]$year),".", unique(eBird_buff_months[[i]]$month))]],
                                               eBird_buff_months[[i]], fun = c("mean"),
                                               progress = FALSE)
  
  eBird_buff_months[[i]]$tmax <- exact_extract(daymet_tmax[[paste0("y", unique(eBird_buff_months[[i]]$year),".", unique(eBird_buff_months[[i]]$month))]],
                                               eBird_buff_months[[i]], fun = c("mean"),
                                               progress = FALSE)
  
  progress(which(names(eBird_buff_months) == i), max.value = length(eBird_buff_months))  
  
}


# Combine to make an object for monthly covariates

monthly_vars <- do.call(rbind, eBird_buff_months)

# Attach and expand to checklists

checklist_vars <- CBCH_eBird %>% 
  select(checklist_id, locality_id, year_lc) %>% 
  inner_join(yearly_vars, by = c("locality_id", "year_lc")) %>% 
  select(-locality_id, -year_lc)

checklist_vars <- CBCH_eBird %>%
  mutate(month = month(observation_date)) %>%
  select(checklist_id, locality_id, year_lc, month) %>% 
  left_join(monthly_vars, by = c("locality_id", "year_lc", "month")) %>% 
  select(-locality_id, -year, -year_lc, -month, -geometry) %>%
  inner_join(checklist_vars, by = "checklist_id")

# Save to CSV, dropping any rows with missing variables

write_csv(drop_na(checklist_vars), 
          "./Data/Cleaned/Environmental-Variables/environmental-variables_checklists.csv", 
          na = "")


## Filling in prediction grid!

# Generate neighborhoods for the prediction grid cell centers

buffers_ps <- as.data.frame(r_3km, cells = TRUE, xy = TRUE) %>% 
  select(cell_id = cell, x, y) %>% 
  st_as_sf(coords = c("x", "y"), crs = laea_crs, remove = FALSE) %>% 
  st_buffer(set_units(1.5, "km")) %>%
  st_transform(crs = 4326)

# Estimate landscape metrics for each cell in the 3km prediction grid

lsm_ps <- NULL

for (i in seq_len(nrow(buffers_ps))) {
  
  buffer_i <- st_transform(buffers_ps[i, ], crs = crs(landcover))
  
  # Crop and mask landcover raster so all values outside buffer are missing
  lsm_ps[[i]] <- crop(landcover[["y2021"]], buffer_i) %>% 
    mask(buffer_i) %>% 
    # calculate landscape metrics
    calculate_lsm(level = "class", metric = c("pland", "ed")) %>% 
    # add variable to uniquely identify each point
    mutate(cell_id = buffer_i$cell_id)
  
  progress(i, nrow(buffers_ps))
  
}

lsm_ps <- bind_rows(lsm_ps) %>% 
  select(cell_id, class, metric, value)

# transform to wide format
lsm_wide_ps <- lsm_ps %>% 
  # fill missing classes with zeros
  complete(cell_id,
           class = lc_classes$class,
           metric = c("ed", "pland"),
           fill = list(value = 0)) %>% 
  # bring in more descriptive names
  inner_join(select(lc_classes, class, label), by = "class") %>% 
  # transform from long to wide format
  pivot_wider(values_from = value,
              names_from = c(class, label, metric),
              names_glue = "{metric}_c{str_pad(class, 2, pad = '0')}_{label}",
              names_sort = TRUE,
              values_fill = 0) %>% 
  arrange(cell_id)

# Create a prediction surface object to fill with values.

pred_surface <- lsm_wide_ps

# Calculate mean and standard deviation of elevation in each grid cell

buffers_ps <- st_transform(buffers_ps, crs = crs(elevation))

pred_surface <- exact_extract(elevation, buffers_ps, fun = c("mean", "stdev"),
                              progress = TRUE) %>% 
  # add variables to uniquely identify each point
  mutate(cell_id = buffers_ps$cell_id) %>% 
  select(cell_id, elevation = mean, elevation_sd = stdev) %>%
  left_join(pred_surface, ., by = "cell_id")

# Calculate other elevation-related values in each grid cell

pred_surface <- exact_extract(slope, buffers_ps, fun = c("mean", "stdev"),
                              progress = TRUE) %>% 
  # add variables to uniquely identify each point
  mutate(cell_id = buffers_ps$cell_id) %>% 
  select(cell_id, slope = mean) %>%
  left_join(pred_surface, ., by = "cell_id")

pred_surface <- exact_extract(eastness, buffers_ps, fun = c("mean", "stdev"),
                              progress = TRUE) %>% 
  # add variables to uniquely identify each point
  mutate(cell_id = buffers_ps$cell_id) %>% 
  select(cell_id, eastness = mean) %>%
  left_join(pred_surface, ., by = "cell_id")

pred_surface <- exact_extract(northness, buffers_ps, fun = c("mean", "stdev"),
                              progress = TRUE) %>% 
  # add variables to uniquely identify each point
  mutate(cell_id = buffers_ps$cell_id) %>% 
  select(cell_id, northness = mean) %>%
  left_join(pred_surface, ., by = "cell_id")

# Calculate mean nighttime light in each grid cell

buffers_ps <- st_transform(buffers_ps, crs = crs(VIIRS))

pred_surface <- exact_extract(VIIRS[["y2021"]], buffers_ps, fun = c("mean"), 
                              progress = TRUE, max_cells_in_memory = 5e+07) %>%
  as.data.frame() %>%
  setNames("mean") %>%
  # add variables to uniquely identify each point
  mutate(cell_id = buffers_ps$cell_id) %>% 
  select(cell_id, nighttime_light = mean) %>%
  left_join(pred_surface, ., by = "cell_id")


# Fetch climate variable values

buffers_ps <- st_transform(buffers_ps, crs = crs(daymet_tmin))

pred_surface <- exact_extract(daymet_prcp[["y2021.5"]], buffers_ps, 
                              fun = c("mean"), progress = TRUE,
                              max_cells_in_memory = 4e+07) %>% 
  as.data.frame() %>%
  setNames("mean") %>%
  # add variables to uniquely identify each point
  mutate(cell_id = buffers_ps$cell_id) %>% 
  select(cell_id, prcp = mean) %>%
  left_join(pred_surface, ., by = "cell_id")


pred_surface <- exact_extract(daymet_tmin[["y2021.5"]], buffers_ps, 
                              fun = c("mean"), progress = TRUE,
                              max_cells_in_memory = 4e+07) %>% 
  as.data.frame() %>%
  setNames("mean") %>%
  # add variables to uniquely identify each point
  mutate(cell_id = buffers_ps$cell_id) %>% 
  select(cell_id, tmin = mean) %>%
  left_join(pred_surface, ., by = "cell_id")

pred_surface <- exact_extract(daymet_tmax[["y2021.5"]], buffers_ps, 
                              fun = c("mean"), progress = TRUE,
                              max_cells_in_memory = 4e+07) %>% 
  as.data.frame() %>%
  setNames("mean") %>%
  # add variables to uniquely identify each point
  mutate(cell_id = buffers_ps$cell_id) %>% 
  select(cell_id, tmax = mean) %>%
  left_join(pred_surface, ., by = "cell_id")

# Fetch road density values

buffers_ps <- st_transform(buffers_ps, crs = st_crs(road_dens))

pred_surface <- exact_extract(road_dens, buffers_ps, fun = c("mean"),
                              progress = TRUE) %>% 
  as.data.frame() %>%
  setNames("mean") %>%
  # add variables to uniquely identify each point
  mutate(cell_id = buffers_ps$cell_id) %>% 
  select(cell_id, road_dens = mean) %>%
  left_join(pred_surface, ., by = "cell_id")

# Fetch canopy height values


buffers_ps <- st_transform(buffers_ps, crs = crs(canopy_height))

pred_surface <- exact_extract(canopy_height, buffers_ps, fun = c("mean"),
                              progress = TRUE) %>%
  as.data.frame() %>%
  setNames("mean") %>%
  # add variables to uniquely identify each point
  mutate(cell_id = buffers_ps$cell_id) %>% 
  select(cell_id, canopy_height = mean) %>%
  left_join(pred_surface, ., by = "cell_id")

# Attach the xy coordinates of the cell centers

pred_surface <- buffers_ps %>% 
  st_drop_geometry() %>% 
  select(cell_id, x, y) %>% 
  inner_join(pred_surface, by = "cell_id")

# Save as csv, dropping any rows with missing variables

write_csv(drop_na(pred_surface),
          "./Data/Cleaned/Prediction-Grids/environmental-variables_prediction-grid.csv", 
          na = "")
