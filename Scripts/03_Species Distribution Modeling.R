################################ INTRO AND PREP ################################

# Script Title: Modeling Chestnut-backed Chickadee Relative Abundance across 
# Bird Conservation Region 5 from bootstrapped datasets inside and outside of the 
# range of the Black-capped Chickadee.

## Install required packages (if not already installed).

librarian::shelf(ebirdst, tidyverse, gridExtra, mccf1, ranger, scam,
                 sf, terra, fields, pdp, raster, precrec, 
                 PresenceAbsence, grid, ggpubr , knitr,
                 mapview, svMisc, doParallel, kableExtra)

# Set number of cores for parallel processing

cores <- detectCores() - 1

# Resolve namespace conflicts

select <- dplyr::select
map <- purrr::map
projection <- raster::projection

# Function to calculate partial dependence for a single predictor

calculate_pd <- function(predictor, er_model, calibration_model,
                         data, x_res = 25, n = 1000) {
  # create prediction grid using quantiles
  x_grid <- quantile(data[[predictor]],
                     probs = seq(from = 0, to = 1, length = x_res),
                     na.rm = TRUE)
  # remove duplicates
  x_grid <- x_grid[!duplicated(signif(x_grid, 8))]
  x_grid <- unname(unique(x_grid))
  grid <- data.frame(predictor = predictor, x = x_grid)
  names(grid) <- c("predictor", predictor)
  
  # subsample training data
  n <- min(n, nrow(data))
  data <- data[sample(seq.int(nrow(data)), size = n, replace = FALSE), ]
  
  # drop focal predictor from data
  data <- data[names(data) != predictor]
  grid <- merge(grid, data, all = TRUE)
  
  # predict encounter rate
  p <- predict(er_model, data = grid)
  
  # summarize
  pd <- grid[, c("predictor", predictor)]
  names(pd) <- c("predictor", "x")
  pd$encounter_rate <- p$predictions[, 2]
  pd <- dplyr::group_by(pd, predictor, x)
  pd <- dplyr::summarise(pd,
                         encounter_rate = mean(encounter_rate, na.rm = TRUE),
                         .groups = "drop")
  
  # calibrate
  pd$encounter_rate <- predict(calibration_model, 
                               newdata = data.frame(pred = pd$encounter_rate), 
                               type = "response")
  pd$encounter_rate <- as.numeric(pd$encounter_rate)
  # constrain to 0-1
  pd$encounter_rate[pd$encounter_rate < 0] <- 0
  pd$encounter_rate[pd$encounter_rate > 1] <- 1
  
  return(pd)
  
}

# Load in habitat variables.

habitat <- read_csv("./Data/Cleaned/Environmental-Variables/environmental-variables_checklists.csv")

# Zero-filled eBird data combined with habitat data

BCCH <- read_csv("./Data/Cleaned/eBird-Basic-Dataset/Zero-Filled-and-Filtered/BCCH/ebd_BCCH_BCR5_zf.csv") %>% 
  inner_join(habitat, by = "checklist_id") %>%
  filter(month(observation_date) %in% c(4,5,6,7), latitude <= 51, is.na(time_observations_started) == FALSE)

CBCH <- read_csv("./Data/Cleaned/eBird-Basic-Dataset/Zero-Filled-and-Filtered/CBCH/ebd_CBCH_BCR5_zf.csv") %>% 
  inner_join(habitat, by = "checklist_id") %>%
  filter(month(observation_date) %in% c(4,5,6,7), latitude <= 51, is.na(time_observations_started) == FALSE)

# Prediction grids

pred_grid <- read_csv("./Data/Cleaned/Prediction-Grids/environmental-variables_prediction-grid.csv")

r <- rast("./Data/Cleaned/Prediction-Grids/prediction-grid_3km.tif")

crs <- crs(r)

# Trim prediction grids to 51 latitude

r <- project(r, "+proj=longlat +datum=WGS84 +no_defs +type=crs") %>% 
  crop(., ext(-152.730904527904, -117.660023882018, 37.3609706347002, 51))

pred_grid <- pred_grid %>% 
  st_as_sf(coords = c("x","y"), crs = crs, remove = FALSE) %>% 
  st_transform(crs = crs(r)) %>% 
  st_crop(., r)

r <- project(r, crs) %>%
  crop(., ext(-2e+05, 6e+05, -1.5e+06, 0))

## Load and project GIS data

# Load in land boundary, project

ne_land <- read_sf("./Data/Raw/Natural-Earth/ne_data.gpkg", "ne_land") %>% 
  st_transform(crs = crs) %>% 
  st_geometry()

# Load in bcr boundary

bcr <- file.path("./Data/Raw/BCR_Terrestrial/BCR_Terrestrial_master_International.shp") %>% 
  read_sf() %>%
  select(bcr_code = BCR, bcr_name = Label) %>% 
  filter(bcr_code == 5) %>%
  st_transform(crs = crs) %>%
  st_geometry()

# Load in country lines, project

ne_country_lines <- read_sf("./Data/Raw/Natural-Earth/ne_data.gpkg", "ne_country_lines") %>% 
  st_transform(crs = crs) %>% 
  st_geometry()

# Load in state/province lines, project

ne_state_lines <- read_sf("./Data/Raw/Natural-Earth/ne_data.gpkg", "ne_state_lines") %>% 
  st_transform(crs = crs) %>% 
  st_geometry()

# Load in detailed coastlines

ne_coast_lines <- read_sf("./Data/Raw/Natural-Earth/ne_10m_coastline/ne_10m_coastline.shp") %>%
  st_transform(crs = crs) %>%
  st_geometry()

# Read in county/regional district lines

ne_counties <- read_sf("./Data/Raw/Natural-Earth/ne_10m_admin_2_counties/ne_10m_admin_2_counties.shp") %>%
  st_transform(crs = crs)

regional_districts <- read_sf("./Data/Raw/BC_Regional-Districts/EBC_REGIONAL_DISTRICTS_SP.geojson") %>%
  st_transform(crs = crs)

# Group counties by our desired subunits.

wa_counties <- c("Snohomish", "King", "Pierce")

or_counties <- c("Multnomah", "Washington", "Clackamas")

mv_counties <- c("Metro Vancouver")

vi_counties <- c("Capital", "Cowichan Valley", "Nanaimo", "Comox Valley", "Alberni-Clayoquot")

sf_counties <- c("Napa", "Marin", "Solano", "Alameda", "San Mateo", "San Francisco", "Contra Costa", "Santa Clara")

nc_counties <- c("Humboldt", "Colusa", "Glenn", "Mendocino", "Trinity", "Tehama", "Shasta", "Del Norte", "Siskiyou", "Sonoma", "Lake")

r_wa <- mask(r, vect(ne_counties[ne_counties$REGION == "WA" & ne_counties$NAME %in% wa_counties,]), touches = TRUE)

r_or <- mask(r, vect(ne_counties[ne_counties$REGION == "OR" & ne_counties$NAME %in% or_counties,]), touches = TRUE)

r_mv <- mask(r, vect(regional_districts[regional_districts$REGIONAL_DISTRICT_NAME == mv_counties,]), touches = TRUE)

r_vi <- mask(r, vect(regional_districts[regional_districts$REGIONAL_DISTRICT_NAME %in% vi_counties,]), touches = TRUE)

r_sf <- mask(r, vect(ne_counties[ne_counties$REGION == "CA" & ne_counties$NAME %in% sf_counties,]), touches = TRUE)

r_nc <- mask(r, vect(ne_counties[ne_counties$REGION == "CA" & ne_counties$NAME %in% nc_counties,]), touches = TRUE)

sf_counties_sp <- ne_counties %>%
  filter(REGION == "CA", NAME %in% sf_counties) %>%
  st_transform(crs = crs) %>%
  st_union()

study_region <- st_union(bcr, sf_counties_sp)

# Read in the Black-capped Chickadee range polygon and trim to BCR-5.

BCCH_range <- st_read("./Data/Raw/eBird-Status-and-Trends/bkcchi_range_2021.gpkg", layer = "range", crs = 4326) %>%
  st_transform(crs = crs) %>%
  st_intersection(., study_region) %>%
  st_buffer(5000)

r_in <- mask(r, BCCH_range, touches = TRUE)

poly_in <- st_as_sf(as.polygons(r_in)) %>%
  st_transform(st_crs(pred_grid))

in_grid <- pred_grid[poly_in,] %>%
  st_drop_geometry()

pred_grid <- st_drop_geometry(pred_grid)

out_grid <- pred_grid %>%
  filter(!(cell_id %in% in_grid$cell_id))


wa_grid <- pred_grid %>% 
  vect(geom = c("x", "y"), crs = crs, keepgeom = TRUE) %>%
  terra::intersect(as.polygons(r_wa), .) %>%
  as.data.frame() %>%
  select(-study_region) %>%
  mutate(county = "WA")

or_grid <- pred_grid %>% 
  vect(geom = c("x", "y"), crs = crs, keepgeom = TRUE) %>%
  terra::intersect(as.polygons(r_or), .) %>%
  as.data.frame() %>%
  select(-study_region) %>%
  mutate(county = "OR")

mv_grid <- pred_grid %>% 
  vect(geom = c("x", "y"), crs = crs, keepgeom = TRUE) %>%
  terra::intersect(as.polygons(r_mv), .) %>%
  as.data.frame() %>%
  select(-study_region) %>%
  mutate(county = "MV")

vi_grid <- pred_grid %>% 
  vect(geom = c("x", "y"), crs = crs, keepgeom = TRUE) %>%
  terra::intersect(as.polygons(r_vi), .) %>%
  as.data.frame() %>%
  select(-study_region) %>%
  mutate(county = "VI")

sf_grid <- pred_grid %>% 
  vect(geom = c("x", "y"), crs = crs, keepgeom = TRUE) %>%
  terra::intersect(as.polygons(r_sf), .) %>%
  as.data.frame() %>%
  select(-study_region) %>%
  mutate(county = "SF")

nc_grid <- pred_grid %>% 
  vect(geom = c("x", "y"), crs = crs, keepgeom = TRUE) %>%
  terra::intersect(as.polygons(r_nc), .) %>%
  as.data.frame() %>%
  select(-study_region) %>%
  mutate(county = "NC")

pred_grid_counties <- rbind(wa_grid, or_grid, mv_grid, vi_grid, sf_grid, nc_grid)

# Read in the Chestnut-backed Chickadee range polygon.

CBCH_range <- st_read("./Data/Raw/eBird-Status-and-Trends/chbchi_range_2021.gpkg", layer = "range", crs = 4326) %>%
  st_transform(crs = crs)

r <- rast(r)

## To complete any of this modeling, we need to select habitat covariates to use.
## We will first fix Median Elevation for inclusion, as well as latitude and 
## longitude, and urban cover. For other habitat covariates, we want to select those with >10%
## cover in our study area, as not to include unnecessary variables. Thus, we
## calculate the percent cover for each of those from our prediction surface.


# Set latitude, elevation, nighttime_light, and cropland for inclusion.

habs_covs <- c("elevation", "northness", "eastness", 
               "slope", "prcp", "tmin", "tmax", "nighttime_light", "road_dens",
               "canopy_height", "pland_c13_urban", "ed_c13_urban", "pland_c12_cropland",
               "ed_c12_cropland")

# Summarize the proportion cover for each habitat variable for each year.

hab_cover <- as.data.frame(matrix(nrow = length(names(select(habitat, starts_with("pland")))),
                                  ncol = 2))

names(hab_cover) <- c("Habitat", "Cover")

hab_cover$Habitat <- names(select(habitat, starts_with("pland")))

for(i in hab_cover$Habitat) {
  
  hab_cover$Cover[hab_cover$Habitat == i] <- sum(pred_grid[,i]/(nrow(pred_grid)*100))*100
  
}

# Remove "pland_" so that we can include the edge density for the desired land 
# cover types as well.

hab_cover$Habitat <- gsub(x = hab_cover$Habitat, pattern = "^pland_", replacement = "")

# Append habitats that pass the 1% cutoff to our habitat covariate object.

habs_covs <- c(habs_covs, names(select(habitat, 
                                       ends_with(hab_cover$Habitat[hab_cover$Cover > 1])))[names(select(habitat, 
                                                                                                        ends_with(hab_cover$Habitat[hab_cover$Cover > 1]))) %in% habs_covs == FALSE])

### Mark CBCH set as in or out of the BCCH range. We will refer to CBCH data 
### outside the range of BCCH as the "out" data, and CBCH data inside the range of 
### BCCH as the "in" data.

CBCH_in_sf <- CBCH %>%
  select(checklist_id, latitude, longitude) %>%
  st_as_sf(crs = 4326, coords = c("longitude", "latitude")) %>%
  st_transform(crs = crs) %>%
  st_intersection(.,BCCH_range)

CBCH$BCCH_range <- ifelse(CBCH$checklist_id %in% CBCH_in_sf$checklist_id,
                          "In", "Out")

CBCH_out_sf <- CBCH %>% 
  filter(BCCH_range == "Out") %>%
  st_as_sf(crs = 4326, coords = c("longitude", "latitude"))

# We want to be able to compare model performance and response both in and out
# of the BCCH range, so we will focus on modeling those covariate values that
# appear in both areas. This also improved model fit when tested by removing some
# extreme outliers.

# Nighttime Light

CBCH <- CBCH[CBCH$nighttime_light <= max(CBCH$nighttime_light[CBCH$BCCH_range == "Out"]),]
BCCH <- BCCH[BCCH$nighttime_light <= max(CBCH$nighttime_light[CBCH$BCCH_range == "Out"]),]

# Road Density

CBCH <- CBCH[CBCH$road_dens <= quantile(CBCH$road_dens[CBCH$BCCH_range == "Out"], 0.95),]
BCCH <- BCCH[BCCH$road_dens <= quantile(CBCH$road_dens[CBCH$BCCH_range == "Out"], 0.95),]

############################## BOOTSTRAPPING LOOP ##############################

# Generate a list of seeds to use for repeatable bootstrapping

seeds <- sample(1:500, size = 100, replace = FALSE)

# If using the seeds used for analysis in Macklin & Jankowski (2025), use this
# code to import our seeds rather than randomly generating them as above.

# seeds <- read_csv("./Data/Cleaned/Seeds/seeds.csv") %>%
#     t() %>%
#     unname() %>%
#     c()

# Initialize a couple of lists to fill.

counties <- c("wa", "or", "mv", "vi", "sf", "nc")

out_assessment <- list()
in_assessment <- list()
all_assessment <- list()
counties_assessment <- list()

out_r_pred <- list()
in_r_pred <- list()
all_r_pred <- list()
counties_r_pred <- list()

for(i in counties) {
  
  counties_assessment[[i]] <- list()
  counties_r_pred[[i]] <- list()
  
}

n <- list()

out_in_tests <- list()

save_path <- "./Scripts/RData/SpeciesDistributionModeling.RData"

# Loop!

for(i in seq_along(seeds)) {
  
  gc()
  
  print(paste0("Working on the ", i, "th bootstrap."))
  
  # Set seed for repeatability
  
  set.seed(seeds[i])
  
  ## Prep eBird Data for Analysis
  
  # Spatiotemporally subset to reduce bias
  
  # Sample one checklist per 3km x 3km x 1 week grid for 
  # each year. Sample detection/non-detection independently.
  
  BCCH_ss <- grid_sample_stratified(BCCH,
                                    res = c(3000, 3000, 7),
                                    obs_column = "species_observed",
                                    sample_by = "type")
  
  CBCH_ss <- grid_sample_stratified(CBCH,
                                    res = c(3000, 3000, 7),
                                    obs_column = "species_observed",
                                    sample_by = "type")
  
  
  # Check that the amount of data set aside for testing is approximately equal
  # for data in and out of the BCCH range.
  
  CBCH_ss$type[CBCH_ss$BCCH_range == "In"] <- if_else(runif(nrow(CBCH_ss[CBCH_ss$BCCH_range == "In",])) <= 0.8, "train", "test")
  
  CBCH_ss$type[CBCH_ss$BCCH_range == "Out"] <- if_else(runif(nrow(CBCH_ss[CBCH_ss$BCCH_range == "Out",])) <= 0.8, "train", "test")
  
  # Build full modeling set to split into independent training and testing tests,
  # with counts, effort covariates, and habitat covariates. 
  
  BCCH_train <- BCCH_ss %>% 
    filter(type == "train") %>% 
    # select only the columns to be used in the model
    select(species_observed, observation_count,
           year, day_of_year, time_observations_started,
           duration_minutes, effort_distance_km,
           number_observers, 
           all_of(habs_covs)) %>%
    filter(is.na(time_observations_started) == FALSE)
  
  CBCH_train <- CBCH_ss %>% 
    filter(type == "train") %>% 
    # select only the columns to be used in the model
    select(species_observed, observation_count,
           year, day_of_year, time_observations_started,
           duration_minutes, effort_distance_km,
           number_observers, 
           all_of(habs_covs), BCCH_range) %>%
    filter(is.na(time_observations_started) == FALSE)
  
  # We'll open this object here to fill as we go.
  
  dataset_sizes <- data.frame(dataset = c("out", "all", "in", counties),
                              train = NA, test = NA)
  
  ############################# MODEL 1: CBCH OUT ###############################
  
  ## First, we will construct model (1): our CBCH model based on data outside of
  ## the BCCH range.
  
  # Filter our training set to those points outside of the BCCH range.
  
  CBCH_out_train <- CBCH_train[CBCH_train$BCCH_range == "Out" & is.na(CBCH_train$time_observations_started) == FALSE,]
  
  # Build step one of our two step hurdle model - modeling encounter rate.
  
  # Calculate detection frequency for the balance random forest
  
  CBCH_out_detfreq <- mean(CBCH_out_train$species_observed)
  
  # Train a random forest model for encounter rate
  
  CBCH_out_train_er <- select(CBCH_out_train, -observation_count, -BCCH_range)
  
  CBCH_out_er_model <- ranger(formula =  as.factor(species_observed) ~ ., 
                              data = CBCH_out_train_er,
                              importance = "permutation",
                              probability = TRUE,
                              replace = TRUE,
                              sample.fraction = c(CBCH_out_detfreq, CBCH_out_detfreq),
                              num.threads = cores)
  
  # Select the mcc-f1 optimizing occurrence threshold
  
  CBCH_out_obs_pred <- tibble(obs = as.integer(CBCH_out_train_er$species_observed), 
                              pred = CBCH_out_er_model$predictions[, 2])
  
  CBCH_out_mcc_f1 <- mccf1(response = CBCH_out_obs_pred$obs,
                           predictor = CBCH_out_obs_pred$pred)
  
  CBCH_out_threshold <- summary(CBCH_out_mcc_f1)$best_threshold[1]
  
  # Calibration model
  
  CBCH_out_calibration_model <- scam(obs ~ s(pred, k = 6, bs = "mpi"), 
                                     gamma = 2,
                                     data = CBCH_out_obs_pred)
  
  # Build step two of our two step hurdle model - modeling expected count.
  
  # Attach the predicted encounter rate based on out of bag samples
  
  CBCH_out_train_count <- CBCH_out_train %>%
    select(-BCCH_range)
  
  CBCH_out_train_count$pred_er <- CBCH_out_er_model$predictions[, 2]
  
  # Subset to only observed or predicted detections
  
  CBCH_out_train_count <- CBCH_out_train_count %>% 
    filter(!is.na(observation_count),
           observation_count > 0 | pred_er > CBCH_out_threshold) %>% 
    select(-species_observed, -pred_er)
  
  # Model counts!
  
  CBCH_out_count_model <- ranger(formula = observation_count ~ .,
                                 data = CBCH_out_train_count,
                                 importance = "permutation",
                                 replace = TRUE,
                                 quantreg = TRUE,
                                 num.threads = cores)
  
  # Check model performance
  
  # Get the test set held out from training
  
  CBCH_out_test <- filter(CBCH_ss, type == "test", BCCH_range == "Out") %>% 
    mutate(species_observed = as.integer(species_observed))
  
  # Predict on test data using random forest model
  
  CBCH_out_testPred_er <- predict(CBCH_out_er_model, data = CBCH_out_test, type = "response")
  
  # Extract probability of detection
  
  CBCH_out_testPred_er <- CBCH_out_testPred_er$predictions[, 2]
  
  # Convert to binary using the threshold
  
  CBCH_out_testPred_binary <- as.integer(CBCH_out_testPred_er > CBCH_out_threshold)
  
  # Calibrate
  
  CBCH_out_testPred_calibrated <- predict(CBCH_out_calibration_model, 
                                          newdata = data.frame(pred = CBCH_out_testPred_er), 
                                          type = "response") %>% 
    as.numeric()
  
  # Add encounter rate to predict count
  
  CBCH_out_test$predicted_er <- CBCH_out_testPred_er
  
  # Estimate count
  
  CBCH_out_testPred_count <- predict(CBCH_out_count_model, data = CBCH_out_test, type = "response")
  
  CBCH_out_testPred_count <- CBCH_out_testPred_count$predictions
  
  # Relative abundance is the product of encounter rate and count
  
  CBCH_out_testPred_abundance <- CBCH_out_testPred_calibrated * CBCH_out_testPred_count
  
  # Combine all estimates together and join to actual observations
  
  CBCH_out_testPred_obs <- data.frame(id = seq_along(CBCH_out_testPred_abundance),
                                      # actual detection/non-detection
                                      obs_detection = as.integer(CBCH_out_test$species_observed),
                                      obs_count = CBCH_out_test$observation_count,
                                      # model estimates
                                      pred_count = CBCH_out_testPred_count,
                                      pred_abundance = CBCH_out_testPred_abundance,
                                      # binary detection/on-detection prediction
                                      pred_binary = CBCH_out_testPred_binary,
                                      # calibrated encounter rate
                                      pred_calibrated = CBCH_out_testPred_calibrated) %>%
    # constrain probabilities to 0-1
    mutate(pred_calibrated = pmin(pmax(pred_calibrated, 0), 1))
  
  # Mean squared error (MSE)
  
  CBCH_out_mse <- mean((CBCH_out_testPred_obs$obs_detection - CBCH_out_testPred_obs$pred_calibrated)^2, na.rm = TRUE)
  
  # spearman correlation, based on in range observations only
  
  CBCH_out_spearman_er <- cor(CBCH_out_testPred_obs$pred_calibrated[CBCH_out_testPred_obs$pred_binary > 0], 
                              CBCH_out_testPred_obs$obs_detection[CBCH_out_testPred_obs$pred_binary > 0], 
                              method = "spearman")
  
  # Precision-recall AUC
  
  CBCH_out_em <- precrec::evalmod(scores = CBCH_out_testPred_obs$pred_binary, 
                                  labels = CBCH_out_testPred_obs$obs_detection)
  
  CBCH_out_pr_auc <- precrec::auc(CBCH_out_em) %>% 
    filter(curvetypes == "PRC") %>% 
    pull(aucs)
  
  CBCH_out_roc_auc <- precrec::auc(CBCH_out_em) %>% 
    filter(curvetypes == "ROC") %>% 
    pull(aucs)
  
  # Calculate metrics for binary prediction: kappa, sensitivity, specificity
  
  CBCH_out_pa_metrics <- CBCH_out_testPred_obs %>% 
    select(id, obs_detection, pred_binary) %>% 
    presence.absence.accuracy(na.rm = TRUE, st.dev = FALSE)
  
  # MCC and F1
  
  CBCH_out_mcc_f1_test <- calculate_mcc_f1(CBCH_out_testPred_obs$obs_detection, 
                                           CBCH_out_testPred_obs$pred_binary)
  
  # subset to only those checklists where detection is predicted
  
  CBCH_out_testPred_detections <- filter(CBCH_out_testPred_obs, pred_binary > 0, 
                                         is.na(obs_count) == FALSE)
  
  # Count metrics
  
  CBCH_out_spearman_count <- cor(CBCH_out_testPred_detections$pred_count, 
                                 CBCH_out_testPred_detections$obs_count,
                                 method = "spearman")
  
  CBCH_out_pearson_log_count <- cor(log(CBCH_out_testPred_detections$pred_count + 1),
                                    log(CBCH_out_testPred_detections$obs_count + 1),
                                    method = "pearson")
  
  # Abundance metrics
  
  CBCH_out_spearman_abundance <- cor(CBCH_out_testPred_detections$pred_abundance, 
                                     CBCH_out_testPred_detections$obs_count,
                                     method = "spearman")
  
  CBCH_out_pearson_log_abundance <- cor(log(CBCH_out_testPred_detections$pred_abundance + 1),
                                        log(CBCH_out_testPred_detections$obs_count + 1),
                                        method = "pearson")
  
  # Combine metrics together
  
  out_assessment[[i]] <- tibble(
    MSE = CBCH_out_mse,
    ER_Spearman = CBCH_out_spearman_er,
    Sensitivity = CBCH_out_pa_metrics$sensitivity,
    Specificity = CBCH_out_pa_metrics$specificity,
    Kappa = CBCH_out_pa_metrics$Kappa,
    PR_AUC = CBCH_out_pr_auc,
    ROC_AUC = CBCH_out_roc_auc,
    MCC = CBCH_out_mcc_f1_test$mcc,
    F1 = CBCH_out_mcc_f1_test$f1,
    Count_Spearman = CBCH_out_spearman_count,
    Log_Count_Pearson = CBCH_out_pearson_log_count,
    Abundance_Spearman = CBCH_out_spearman_abundance,
    Log_Abundance_Pearson = CBCH_out_pearson_log_abundance
  )
  
  # Now, predict!
  
  # Find peak time of day for CBCH detection from partial dependence
  
  CBCH_out_pd_time <- calculate_pd("time_observations_started",
                                   er_model = CBCH_out_er_model, 
                                   data = CBCH_out_train,
                                   calibration_model = CBCH_out_calibration_model,
                                   # make estimates at 30 minute intervals
                                   # using a subset of the training dataset
                                   x_res = 2 * 24, n = 1000) %>% 
    select(time_observations_started = x, encounter_rate)
  
  # hours with at least 1% of checklists
  
  CBCH_out_search_times <- CBCH_out_train %>% 
    mutate(time_observations_started = floor(time_observations_started)) %>%
    count(time_observations_started) %>% 
    mutate(pct = n / sum(n)) %>% 
    filter(pct >= 0.01)
  
  CBCH_out_t_peak <- CBCH_out_pd_time %>% 
    filter(floor(time_observations_started) %in% CBCH_out_search_times$time_observations_started) %>% 
    slice_max(order_by = encounter_rate) %>% 
    pull(time_observations_started)
  
  # Add effort covariates to prediction grid
  
  CBCH_out_pred_grid_eff <- pred_grid %>% 
    mutate(observation_date = ymd("2022-06-01"),
           year = year(observation_date),
           day_of_year = yday(observation_date),
           time_observations_started = CBCH_out_t_peak,
           duration_minutes = 120,
           effort_distance_km = 1,
           number_observers = 1)
  
  
  # Predict encounter rate first.
  
  CBCH_out_pred_er <- predict(CBCH_out_er_model, data = CBCH_out_pred_grid_eff, 
                              type = "response")
  
  CBCH_out_pred_er <- as.data.frame(cbind(CBCH_out_pred_grid_eff$cell_id,
                                          CBCH_out_pred_er$predictions[,2])) %>%
    setNames(c("cell_id", "predicted_er"))
  
  # Apply threshold
  
  CBCH_out_pred_er$predicted_binary <- as.integer(CBCH_out_pred_er$predicted_er > CBCH_out_threshold)
  
  # Apply calibration
  
  CBCH_out_pred_er$calibrated_er <- predict(CBCH_out_calibration_model, 
                                            data.frame(pred = CBCH_out_pred_er$predicted_er), 
                                            type = "response") %>% 
    as.numeric()
  
  # Add predicted encounter rate required for count estimates
  
  CBCH_out_pred_grid_eff <- left_join(CBCH_out_pred_grid_eff, 
                                      CBCH_out_pred_er, 
                                      by = "cell_id")
  
  # Predict count estimate
  
  CBCH_out_pred_count <- predict(CBCH_out_count_model, 
                                 data = CBCH_out_pred_grid_eff, 
                                 type = "response")
  
  CBCH_out_pred_count <- as.data.frame(cbind(CBCH_out_pred_grid_eff$cell_id,
                                             CBCH_out_pred_count$predictions)) %>%
    setNames(c("cell_id", "count"))
  
  
  # Add estimates to prediction grid
  
  CBCH_out_pred_grid_eff <- left_join(CBCH_out_pred_grid_eff, 
                                      CBCH_out_pred_count,
                                      by = "cell_id")
  
  # Bind predictions together
  
  CBCH_out_predictions <- CBCH_out_pred_grid_eff %>% 
    mutate(in_range = predicted_binary, encounter_rate = calibrated_er) %>%
    select(cell_id, x, y, in_range, encounter_rate, count) %>% 
    mutate(encounter_rate = pmin(pmax(encounter_rate, 0), 1))
  
  # Calculate relative abundance and rasterize.
  
  out_r_pred[[i]] <- CBCH_out_predictions %>% 
    
    # Estimate relative abundance
    mutate(abundance = encounter_rate * count) %>% 
    
    # Convert to spatial features
    st_as_sf(coords = c("x", "y"), crs = crs) %>% 
    select(cell_id, in_range, encounter_rate, count, abundance) %>%
    
    # Rasterize
    rasterize(r, field = c("cell_id", "in_range", "encounter_rate", "count", "abundance"),
              fun = "mean") %>% 
    setNames(c("cell_id", "in_range", "encounter_rate", "count", "abundance")) %>%
    as.data.frame()

  # Store dataset sizes
  
  dataset_sizes$train[dataset_sizes$dataset == "out"] <- nrow(CBCH_out_train)
  dataset_sizes$test[dataset_sizes$dataset == "out"] <- nrow(CBCH_out_test)
  
  print(paste0("Bootstrap ", i, ": Out Model Done (1/", 3+length(counties), ")"))
  
  ############################# MODEL 2: CBCH IN #################################
  
  ## Next, we will construct model (2): our CBCH model based on data inside of
  ## the BCCH range.
  
  # Filter our training set to those points inside of the BCCH range.
  
  CBCH_in_train <- CBCH_train[CBCH_train$BCCH_range == "In" & is.na(CBCH_train$time_observations_started) == FALSE,]
  
  # First, we need to define a BCCH range.
  
  # Build step one of our two step hurdle model for BCCH - modeling encounter rate.
  
  # Calculate detection frequency for the balance random forest
  
  BCCH_detfreq <- mean(BCCH_train$species_observed)
  
  # Train a random forest model for encounter rate
  
  BCCH_train_er <- select(BCCH_train, -observation_count) %>%
    filter(is.na(time_observations_started) == FALSE)
  
  BCCH_er_model <- ranger(formula =  as.factor(species_observed) ~ ., 
                          data = BCCH_train_er,
                          importance = "permutation",
                          probability = TRUE,
                          replace = TRUE,
                          sample.fraction = c(BCCH_detfreq, BCCH_detfreq),
                          num.threads = cores)
  
  # Select the mcc-f1 optimizing occurrence threshold
  
  BCCH_obs_pred <- tibble(obs = as.integer(BCCH_train_er$species_observed), 
                          pred = BCCH_er_model$predictions[, 2])
  
  BCCH_mcc_f1 <- mccf1(response = BCCH_obs_pred$obs, predictor = BCCH_obs_pred$pred)
  
  BCCH_threshold <- summary(BCCH_mcc_f1)$best_threshold[1]
  
  # Now, predict from our first model the BCCH abundance at each checklist in the
  # range of BCCH.
  
  # First, predict encounter rate
  
  BCCH_pred_er <- predict(BCCH_er_model,
                          CBCH_in_train,
                          type = "response")
  
  BCCH_pred_er <- as.data.frame(BCCH_pred_er$predictions[,2]) %>%
    setNames("predicted_er")
  
  # Apply threshold
  
  BCCH_pred_er$predicted_binary <- as.integer(BCCH_pred_er$predicted_er > BCCH_threshold)
  
  # Join to training data
  
  CBCH_in_train$BCCH_range_binary <- BCCH_pred_er$predicted_binary
  
  # Now that we have estimates of BCCH range at each point, build model (2)
  
  # Filter training data to cells with BCCH predicted to occur
  
  CBCH_in_polygon_train <- CBCH_in_train
  
  CBCH_in_train <- CBCH_in_train %>%
    filter(BCCH_range_binary == 1) %>%
    select(-BCCH_range, -BCCH_range_binary)
  
  # Build step one of our two step hurdle model - modeling encounter rate.
  
  # Calculate detection frequency for the balance random forest
  
  CBCH_in_detfreq <- mean(CBCH_in_train$species_observed)
  
  # Train a random forest model for encounter rate
  
  CBCH_in_train_er <- select(CBCH_in_train, -observation_count)
  
  CBCH_in_er_model <- ranger(formula =  as.factor(species_observed) ~ ., 
                             data = CBCH_in_train_er,
                             importance = "permutation",
                             probability = TRUE,
                             replace = TRUE,
                             sample.fraction = c(CBCH_in_detfreq, CBCH_in_detfreq),
                             num.threads = cores)
  
  # Select the mcc-f1 optimizing occurrence threshold
  
  CBCH_in_obs_pred <- tibble(obs = as.integer(CBCH_in_train_er$species_observed), 
                             pred = CBCH_in_er_model$predictions[, 2])
  
  CBCH_in_mcc_f1 <- mccf1(response = CBCH_in_obs_pred$obs,
                          predictor = CBCH_in_obs_pred$pred)
  
  CBCH_in_threshold <- summary(CBCH_in_mcc_f1)$best_threshold[1]
  
  # Calibration model
  
  CBCH_in_calibration_model <- scam(obs ~ s(pred, k = 6, bs = "mpi"), 
                                    gamma = 2,
                                    data = CBCH_in_obs_pred)
  
  # Build step two of our two step hurdle model - modeling expected count.
  
  # Attach the predicted encounter rate based on in of bag samples
  
  CBCH_in_train_count <- CBCH_in_train
  
  CBCH_in_train_count$pred_er <- CBCH_in_er_model$predictions[, 2]
  
  # Subset to only observed or predicted detections
  
  CBCH_in_train_count <- CBCH_in_train_count %>% 
    filter(!is.na(observation_count),
           observation_count > 0 | pred_er > CBCH_in_threshold) %>% 
    select(-species_observed, -pred_er)
  
  # Model counts!
  
  CBCH_in_count_model <- ranger(formula = observation_count ~ .,
                                data = CBCH_in_train_count,
                                importance = "permutation",
                                replace = TRUE,
                                quantreg = TRUE, num.threads = cores)
  
  # Check model performance
  
  # Get the test set held in from training
  
  CBCH_in_test <- filter(CBCH_ss, type == "test", BCCH_range == "In") %>% 
    mutate(species_observed = as.integer(species_observed))
  
  # Now, predict from our first model the BCCH binary occurence at each checklist 
  # in the range of BCCH.
  
  # First, predict encounter rate
  
  BCCH_all_pred_er <- predict(BCCH_er_model,
                              CBCH_in_test,
                              type = "response")
  
  BCCH_all_pred_er <- as.data.frame(BCCH_all_pred_er$predictions[,2]) %>%
    setNames("predicted_er")
  
  # Apply threshold
  
  BCCH_all_pred_er$predicted_binary <- as.integer(BCCH_all_pred_er$predicted_er > BCCH_threshold)
  
  # Join to training data
  
  CBCH_in_test$BCCH_range_binary <- BCCH_all_pred_er$predicted_binary
  
  # Filter training data to cells with BCCH predicted to occur
  
  CBCH_in_polygon_test <- CBCH_in_test
  
  CBCH_in_test <- CBCH_in_test %>%
    filter(BCCH_range_binary == 1) %>%
    select(-BCCH_range, -BCCH_range_binary)
  
  # Predict on test data using random forest model
  
  CBCH_in_testPred_er <- predict(CBCH_in_er_model, data = CBCH_in_test, type = "response")
  
  # Extract probability of detection
  
  CBCH_in_testPred_er <- CBCH_in_testPred_er$predictions[, 2]
  
  # Convert to binary using the threshold
  
  CBCH_in_testPred_binary <- as.integer(CBCH_in_testPred_er > CBCH_in_threshold)
  
  # Calibrate
  
  CBCH_in_testPred_calibrated <- predict(CBCH_in_calibration_model, 
                                         newdata = data.frame(pred = CBCH_in_testPred_er), 
                                         type = "response") %>% 
    as.numeric()
  
  # Add encounter rate to predict count
  
  CBCH_in_test$predicted_er <- CBCH_in_testPred_er
  
  # Estimate count
  
  CBCH_in_testPred_count <- predict(CBCH_in_count_model, data = CBCH_in_test, type = "response")
  
  CBCH_in_testPred_count <- CBCH_in_testPred_count$predictions
  
  # Relative abundance is the product of encounter rate and count
  
  CBCH_in_testPred_abundance <- CBCH_in_testPred_calibrated * CBCH_in_testPred_count
  
  # Combine all estimates together and join to actual observations
  
  CBCH_in_testPred_obs <- data.frame(id = seq_along(CBCH_in_testPred_abundance),
                                     # actual detection/non-detection
                                     obs_detection = as.integer(CBCH_in_test$species_observed),
                                     obs_count = CBCH_in_test$observation_count,
                                     # model estimates
                                     pred_count = CBCH_in_testPred_count,
                                     pred_abundance = CBCH_in_testPred_abundance,
                                     # binary detection/on-detection prediction
                                     pred_binary = CBCH_in_testPred_binary,
                                     # calibrated encounter rate
                                     pred_calibrated = CBCH_in_testPred_calibrated) %>%
    # constrain probabilities to 0-1
    mutate(pred_calibrated = pmin(pmax(pred_calibrated, 0), 1))
  
  # Mean squared error (MSE)
  
  CBCH_in_mse <- mean((CBCH_in_testPred_obs$obs_detection - CBCH_in_testPred_obs$pred_calibrated)^2, na.rm = TRUE)
  
  # spearman correlation, based on in range observations only
  
  CBCH_in_spearman_er <- cor(CBCH_in_testPred_obs$pred_calibrated[CBCH_in_testPred_obs$pred_binary > 0], 
                             CBCH_in_testPred_obs$obs_detection[CBCH_in_testPred_obs$pred_binary > 0], 
                             method = "spearman")
  
  # Precision-recall AUC
  
  CBCH_in_em <- precrec::evalmod(scores = CBCH_in_testPred_obs$pred_binary, 
                                 labels = CBCH_in_testPred_obs$obs_detection)
  
  CBCH_in_pr_auc <- precrec::auc(CBCH_in_em) %>% 
    filter(curvetypes == "PRC") %>% 
    pull(aucs)
  
  CBCH_in_roc_auc <- precrec::auc(CBCH_in_em) %>% 
    filter(curvetypes == "ROC") %>% 
    pull(aucs)
  
  # Calculate metrics for binary prediction: kappa, sensitivity, specificity
  
  CBCH_in_pa_metrics <- CBCH_in_testPred_obs %>% 
    select(id, obs_detection, pred_binary) %>% 
    presence.absence.accuracy(na.rm = TRUE, st.dev = FALSE)
  
  # MCC and F1
  
  CBCH_in_mcc_f1_test <- calculate_mcc_f1(CBCH_in_testPred_obs$obs_detection, 
                                          CBCH_in_testPred_obs$pred_binary)
  
  # subset to only those checklists where detection is predicted
  
  CBCH_in_testPred_detections <- filter(CBCH_in_testPred_obs, pred_binary > 0, 
                                        is.na(obs_count) == FALSE)
  
  # Count metrics
  
  CBCH_in_spearman_count <- cor(CBCH_in_testPred_detections$pred_count, 
                                CBCH_in_testPred_detections$obs_count,
                                method = "spearman")
  
  CBCH_in_pearson_log_count <- cor(log(CBCH_in_testPred_detections$pred_count + 1),
                                   log(CBCH_in_testPred_detections$obs_count + 1),
                                   method = "pearson")
  
  # Abundance metrics
  
  CBCH_in_spearman_abundance <- cor(CBCH_in_testPred_detections$pred_abundance, 
                                    CBCH_in_testPred_detections$obs_count,
                                    method = "spearman")
  
  CBCH_in_pearson_log_abundance <- cor(log(CBCH_in_testPred_detections$pred_abundance + 1),
                                       log(CBCH_in_testPred_detections$obs_count + 1),
                                       method = "pearson")
  
  # Combine metrics together
  
  in_assessment[[i]] <- tibble(
    MSE = CBCH_in_mse,
    ER_Spearman = CBCH_in_spearman_er,
    Sensitivity = CBCH_in_pa_metrics$sensitivity,
    Specificity = CBCH_in_pa_metrics$specificity,
    Kappa = CBCH_in_pa_metrics$Kappa,
    PR_AUC = CBCH_in_pr_auc,
    ROC_AUC = CBCH_in_roc_auc,
    MCC = CBCH_in_mcc_f1_test$mcc,
    F1 = CBCH_in_mcc_f1_test$f1,
    Count_Spearman = CBCH_in_spearman_count,
    Log_Count_Pearson = CBCH_in_pearson_log_count,
    Abundance_Spearman = CBCH_in_spearman_abundance,
    Log_Abundance_Pearson = CBCH_in_pearson_log_abundance
  )
  
  # Now, predict!
  
  # Find peak time of day for CBCH detection from partial dependence
  
  CBCH_in_pd_time <- calculate_pd("time_observations_started",
                                  er_model = CBCH_in_er_model,
                                  calibration_model = CBCH_in_calibration_model,
                                  data = CBCH_in_train,
                                  # make estimates at 30 minute intervals
                                  # using a subset of the training dataset
                                  x_res = 2 * 24, n = 1000) %>% 
    select(time_observations_started = x, encounter_rate)
  
  # hours with at least 1% of checklists
  
  CBCH_in_search_times <- CBCH_in_train %>% 
    mutate(time_observations_started = floor(time_observations_started)) %>%
    count(time_observations_started) %>% 
    mutate(pct = n / sum(n)) %>% 
    filter(pct >= 0.01)
  
  CBCH_in_t_peak <- CBCH_in_pd_time %>% 
    filter(floor(time_observations_started) %in% CBCH_in_search_times$time_observations_started) %>% 
    slice_max(order_by = encounter_rate) %>% 
    pull(time_observations_started)
  
  # Add effort covariates to prediction grid
  
  CBCH_in_pred_grid_eff <- pred_grid %>% 
    mutate(observation_date = ymd("2022-06-01"),
           year = year(observation_date),
           day_of_year = yday(observation_date),
           time_observations_started = CBCH_in_t_peak,
           duration_minutes = 120,
           effort_distance_km = 1,
           number_observers = 1)
  
  
  # Predict encounter rate first.
  
  CBCH_in_pred_er <- predict(CBCH_in_er_model, data = CBCH_in_pred_grid_eff, 
                             type = "response")
  
  CBCH_in_pred_er <- as.data.frame(cbind(CBCH_in_pred_grid_eff$cell_id,
                                         CBCH_in_pred_er$predictions[,2])) %>%
    setNames(c("cell_id", "predicted_er"))
  
  # Apply threshold
  
  CBCH_in_pred_er$predicted_binary <- as.integer(CBCH_in_pred_er$predicted_er > CBCH_in_threshold)
  
  # Apply calibration
  
  CBCH_in_pred_er$calibrated_er <- predict(CBCH_in_calibration_model, 
                                           data.frame(pred = CBCH_in_pred_er$predicted_er), 
                                           type = "response") %>% 
    as.numeric()
  
  # Add predicted encounter rate required for count estimates
  
  CBCH_in_pred_grid_eff <- left_join(CBCH_in_pred_grid_eff, 
                                     CBCH_in_pred_er, 
                                     by = "cell_id")
  
  # Predict count estimate
  
  CBCH_in_pred_count <- predict(CBCH_in_count_model, 
                                data = CBCH_in_pred_grid_eff, 
                                type = "response")
  
  CBCH_in_pred_count <- as.data.frame(cbind(CBCH_in_pred_grid_eff$cell_id,
                                            CBCH_in_pred_count$predictions)) %>%
    setNames(c("cell_id", "count"))
  
  
  # Add estimates to prediction grid
  
  CBCH_in_pred_grid_eff <- left_join(CBCH_in_pred_grid_eff, 
                                     CBCH_in_pred_count,
                                     by = "cell_id")
  
  # Bind predictions together
  
  CBCH_in_predictions <- CBCH_in_pred_grid_eff %>% 
    mutate(in_range = predicted_binary, encounter_rate = calibrated_er) %>%
    select(cell_id, x, y, in_range, encounter_rate, count) %>% 
    mutate(encounter_rate = pmin(pmax(encounter_rate, 0), 1))
  
  # Calculate relative abundance and rasterize.
  
  in_r_pred[[i]] <- CBCH_in_predictions %>% 
    
    # Estimate relative abundance
    mutate(abundance = encounter_rate * count) %>% 
    
    # Convert to spatial features
    st_as_sf(coords = c("x", "y"), crs = crs) %>% 
    select(cell_id, in_range, encounter_rate, count, abundance) %>%
    
    # Rasterize
    rasterize(r, field = c("cell_id", "in_range", "encounter_rate", "count", "abundance"),
              fun = "mean") %>% 
    setNames(c("cell_id", "in_range", "encounter_rate", "count", "abundance")) %>%
    as.data.frame()

  # Store dataset sizes
  
  dataset_sizes$train[dataset_sizes$dataset == "in"] <- nrow(CBCH_in_train)
  dataset_sizes$test[dataset_sizes$dataset == "in"] <- nrow(CBCH_in_test)
  
  print(paste0("Bootstrap ", i, ": In Model Done (2/", 3+length(counties), ")"))
  
  ############################ MODEL 3: ALL ###################################
  
  # Combine the two!
  
  CBCH_all_train <- CBCH_out_train %>%
    select(-BCCH_range)
  
  CBCH_all_train <- CBCH_in_train %>%
    rbind(CBCH_all_train)
  
  ### Start building our models!
  
  # Build step one of our two step hurdle model - modeling encounter rate.
  
  # Calculate detection frequency for the balance random forest
  
  CBCH_all_detfreq <- mean(CBCH_all_train$species_observed)
  
  # Train a random forest model for encounter rate
  
  CBCH_all_train_er <- select(CBCH_all_train, -observation_count)
  
  CBCH_all_er_model <- ranger(formula =  as.factor(species_observed) ~ ., 
                              data = CBCH_all_train_er,
                              importance = "permutation",
                              probability = TRUE,
                              replace = TRUE,
                              sample.fraction = c(CBCH_all_detfreq, CBCH_all_detfreq))
  
  # Select the mcc-f1 optimizing occurrence threshold
  
  CBCH_all_obs_pred <- tibble(obs = as.integer(CBCH_all_train_er$species_observed), 
                              pred = CBCH_all_er_model$predictions[, 2])
  
  CBCH_all_mcc_f1 <- mccf1(response = CBCH_all_obs_pred$obs,
                           predictor = CBCH_all_obs_pred$pred)
  
  CBCH_all_threshold <- summary(CBCH_all_mcc_f1)$best_threshold[1]
  
  # Calibration model
  
  CBCH_all_calibration_model <- scam(obs ~ s(pred, k = 6, bs = "mpi"), 
                                     gamma = 2,
                                     data = CBCH_all_obs_pred)
  
  # Build step two of our two step hurdle model - modeling expected count.
  
  # Attach the predicted encounter rate based on in of bag samples
  
  CBCH_all_train_count <- CBCH_all_train
  
  CBCH_all_train_count$pred_er <- CBCH_all_er_model$predictions[, 2]
  
  # Subset to only observed or predicted detections
  
  CBCH_all_train_count <- CBCH_all_train_count %>% 
    filter(!is.na(observation_count),
           observation_count > 0 | pred_er > CBCH_all_threshold) %>% 
    select(-species_observed, -pred_er)
  
  # Model counts!
  
  CBCH_all_count_model <- ranger(formula = observation_count ~ .,
                                 data = CBCH_all_train_count,
                                 importance = "permutation",
                                 replace = TRUE,
                                 quantreg = TRUE, num.threads = cores)
  
  
  ## Create testing set
  
  # Get the test set held out from training for the Out Set
  
  CBCH_out_test <- filter(CBCH_ss, type == "test", BCCH_range == "Out") %>% 
    mutate(species_observed = as.integer(species_observed))
  
  
  # Get the test set held in from training
  
  CBCH_in_test <- filter(CBCH_ss, type == "test", BCCH_range == "In") %>% 
    mutate(species_observed = as.integer(species_observed))
  
  # Now, predict from our first model the BCCH binary occurence at each checklist 
  # in the range of BCCH.
  
  # First, predict encounter rate
  
  BCCH_all_pred_er <- predict(BCCH_er_model,
                              CBCH_in_test,
                              type = "response")
  
  BCCH_all_pred_er <- as.data.frame(BCCH_all_pred_er$predictions[,2]) %>%
    setNames("predicted_er")
  
  # Apply threshold
  
  BCCH_all_pred_er$predicted_binary <- as.integer(BCCH_all_pred_er$predicted_er > BCCH_threshold)
  
  # Join to training data
  
  CBCH_in_test$BCCH_range_binary <- BCCH_all_pred_er$predicted_binary
  
  # Filter training data to cells with BCCH predicted to occur
  
  CBCH_in_polygon_test <- CBCH_in_test
  
  CBCH_in_test <- CBCH_in_test %>%
    filter(BCCH_range_binary == 1) %>%
    select(-BCCH_range, -BCCH_range_binary)
  
  # Make our all test set
  
  CBCH_all_test <- CBCH_out_test %>%
    mutate(BCCH_range = "Out")
  
  CBCH_all_test <- CBCH_in_test %>%
    mutate(BCCH_range = "In") %>%
    rbind(CBCH_all_test)
  
  CBCH_all_test$BCCH_range <- as.factor(CBCH_all_test$BCCH_range)
  
  
  CBCH_all_testPred_er <- predict(CBCH_all_er_model, data = CBCH_all_test, type = "response")
  
  # Extract probability of detection
  
  CBCH_all_testPred_er <- CBCH_all_testPred_er$predictions[, 2]
  
  # Convert to binary using the threshold
  
  CBCH_all_testPred_binary <- as.integer(CBCH_all_testPred_er > CBCH_all_threshold)
  
  # Calibrate
  
  CBCH_all_testPred_calibrated <- predict(CBCH_all_calibration_model, 
                                          newdata = data.frame(pred = CBCH_all_testPred_er), 
                                          type = "response") %>% 
    as.numeric()
  
  # Add encounter rate to predict count
  
  CBCH_all_test$predicted_er <- CBCH_all_testPred_er
  
  # Estimate count
  
  CBCH_all_testPred_count <- predict(CBCH_all_count_model, data = CBCH_all_test, type = "response")
  
  CBCH_all_testPred_count <- CBCH_all_testPred_count$predictions
  
  # Relative abundance is the product of encounter rate and count
  
  CBCH_all_testPred_abundance <- CBCH_all_testPred_calibrated * CBCH_all_testPred_count
  
  # Combine all estimates together and join to actual observations
  
  CBCH_all_testPred_obs <- data.frame(id = seq_along(CBCH_all_testPred_abundance),
                                      # actual detection/non-detection
                                      obs_detection = as.integer(CBCH_all_test$species_observed),
                                      obs_count = CBCH_all_test$observation_count,
                                      # model estimates
                                      pred_count = CBCH_all_testPred_count,
                                      pred_abundance = CBCH_all_testPred_abundance,
                                      # binary detection/on-detection prediction
                                      pred_binary = CBCH_all_testPred_binary,
                                      # calibrated encounter rate
                                      pred_calibrated = CBCH_all_testPred_calibrated) %>%
    # constrain probabilities to 0-1
    mutate(pred_calibrated = pmin(pmax(pred_calibrated, 0), 1))
  
  # Mean squared error (MSE)
  
  CBCH_all_mse <- mean((CBCH_all_testPred_obs$obs_detection - CBCH_all_testPred_obs$pred_calibrated)^2, na.rm = TRUE)
  
  # spearman correlation, based on in range observations only
  
  CBCH_all_spearman_er <- cor(CBCH_all_testPred_obs$pred_calibrated[CBCH_all_testPred_obs$pred_binary > 0], 
                              CBCH_all_testPred_obs$obs_detection[CBCH_all_testPred_obs$pred_binary > 0], 
                              method = "spearman")
  
  # Precision-recall AUC
  
  CBCH_all_em <- precrec::evalmod(scores = CBCH_all_testPred_obs$pred_binary, 
                                  labels = CBCH_all_testPred_obs$obs_detection)
  
  CBCH_all_pr_auc <- precrec::auc(CBCH_all_em) %>% 
    filter(curvetypes == "PRC") %>% 
    pull(aucs)
  
  CBCH_all_roc_auc <- precrec::auc(CBCH_all_em) %>% 
    filter(curvetypes == "ROC") %>% 
    pull(aucs)
  
  # Calculate metrics for binary prediction: kappa, sensitivity, specificity
  
  CBCH_all_pa_metrics <- CBCH_all_testPred_obs %>% 
    select(id, obs_detection, pred_binary) %>% 
    presence.absence.accuracy(na.rm = TRUE, st.dev = FALSE)
  
  # MCC and F1
  
  CBCH_all_mcc_f1_test <- calculate_mcc_f1(CBCH_all_testPred_obs$obs_detection, 
                                           CBCH_all_testPred_obs$pred_binary)
  
  # subset to only those checklists where detection is predicted
  
  CBCH_all_testPred_detections <- filter(CBCH_all_testPred_obs, pred_binary > 0, 
                                         is.na(obs_count) == FALSE)
  
  # Count metrics
  
  CBCH_all_spearman_count <- cor(CBCH_all_testPred_detections$pred_count, 
                                 CBCH_all_testPred_detections$obs_count,
                                 method = "spearman")
  
  CBCH_all_pearson_log_count <- cor(log(CBCH_all_testPred_detections$pred_count + 1),
                                    log(CBCH_all_testPred_detections$obs_count + 1),
                                    method = "pearson")
  
  # Abundance metrics
  
  CBCH_all_spearman_abundance <- cor(CBCH_all_testPred_detections$pred_abundance, 
                                     CBCH_all_testPred_detections$obs_count,
                                     method = "spearman")
  
  CBCH_all_pearson_log_abundance <- cor(log(CBCH_all_testPred_detections$pred_abundance + 1),
                                        log(CBCH_all_testPred_detections$obs_count + 1),
                                        method = "pearson")
  
  # Combine metrics together
  
  all_assessment[[i]] <- tibble(
    MSE = CBCH_all_mse,
    ER_Spearman = CBCH_all_spearman_er,
    Sensitivity = CBCH_all_pa_metrics$sensitivity,
    Specificity = CBCH_all_pa_metrics$specificity,
    Kappa = CBCH_all_pa_metrics$Kappa,
    PR_AUC = CBCH_all_pr_auc,
    ROC_AUC = CBCH_all_roc_auc,
    MCC = CBCH_all_mcc_f1_test$mcc,
    F1 = CBCH_all_mcc_f1_test$f1,
    Count_Spearman = CBCH_all_spearman_count,
    Log_Count_Pearson = CBCH_all_pearson_log_count,
    Abundance_Spearman = CBCH_all_spearman_abundance,
    Log_Abundance_Pearson = CBCH_all_pearson_log_abundance
  )
  
  # Now, predict!
  
  # Find peak time of day for CBCH detection from partial dependence
  
  CBCH_all_pd_time <- calculate_pd("time_observations_started",
                                   er_model = CBCH_all_er_model,
                                   calibration_model = CBCH_all_calibration_model,
                                   data = CBCH_all_train,
                                   # make estimates at 30 minute intervals
                                   # using a subset of the training dataset
                                   x_res = 2 * 24, n = 1000) %>% 
    select(time_observations_started = x, encounter_rate)
  
  # hours with at least 1% of checklists
  
  CBCH_all_search_times <- CBCH_all_train %>% 
    mutate(time_observations_started = floor(time_observations_started)) %>%
    count(time_observations_started) %>% 
    mutate(pct = n / sum(n)) %>% 
    filter(pct >= 0.01)
  
  CBCH_all_t_peak <- CBCH_all_pd_time %>% 
    filter(floor(time_observations_started) %in% CBCH_all_search_times$time_observations_started) %>% 
    slice_max(order_by = encounter_rate) %>% 
    pull(time_observations_started)
  
  # Add effort covariates to prediction grid
  
  CBCH_all_pred_grid_eff <- pred_grid %>% 
    mutate(observation_date = ymd("2022-06-01"),
           year = year(observation_date),
           day_of_year = yday(observation_date),
           time_observations_started = CBCH_all_t_peak,
           duration_minutes = 120,
           effort_distance_km = 1,
           number_observers = 1)
  
  
  # Predict encounter rate first.
  
  CBCH_all_pred_er <- predict(CBCH_all_er_model, data = CBCH_all_pred_grid_eff, 
                              type = "response")
  
  CBCH_all_pred_er <- as.data.frame(cbind(CBCH_all_pred_grid_eff$cell_id,
                                          CBCH_all_pred_er$predictions[,2])) %>%
    setNames(c("cell_id", "predicted_er"))
  
  # Apply threshold
  
  CBCH_all_pred_er$predicted_binary <- as.integer(CBCH_all_pred_er$predicted_er > CBCH_all_threshold)
  
  # Apply calibration
  
  CBCH_all_pred_er$calibrated_er <- predict(CBCH_all_calibration_model, 
                                            data.frame(pred = CBCH_all_pred_er$predicted_er), 
                                            type = "response") %>% 
    as.numeric()
  
  # Add predicted encounter rate required for count estimates
  
  CBCH_all_pred_grid_eff <- left_join(CBCH_all_pred_grid_eff, 
                                      CBCH_all_pred_er, 
                                      by = "cell_id")
  
  # Predict count estimate
  
  CBCH_all_pred_count <- predict(CBCH_all_count_model, 
                                 data = CBCH_all_pred_grid_eff, 
                                 type = "response")
  
  CBCH_all_pred_count <- as.data.frame(cbind(CBCH_all_pred_grid_eff$cell_id,
                                             CBCH_all_pred_count$predictions)) %>%
    setNames(c("cell_id", "count"))
  
  
  # Add estimates to prediction grid
  
  CBCH_all_pred_grid_eff <- left_join(CBCH_all_pred_grid_eff, 
                                      CBCH_all_pred_count,
                                      by = "cell_id")
  
  # Bind predictions together
  
  CBCH_all_predictions <- CBCH_all_pred_grid_eff %>% 
    mutate(in_range = predicted_binary, encounter_rate = calibrated_er) %>%
    select(cell_id, x, y, in_range, encounter_rate, count) %>% 
    mutate(encounter_rate = pmin(pmax(encounter_rate, 0), 1))
  
  # Calculate relative abundance and rasterize.
  
  all_r_pred[[i]] <- CBCH_all_predictions %>% 
    
    # Estimate relative abundance
    mutate(abundance = encounter_rate * count) %>% 
    
    # Convert to spatial features
    st_as_sf(coords = c("x", "y"), crs = crs) %>% 
    select(cell_id, in_range, encounter_rate, count, abundance) %>%
    
    # Rasterize
    rasterize(r, field = c("cell_id", "in_range", "encounter_rate", "count", "abundance"),
              fun = "mean") %>% 
    setNames(c("cell_id", "in_range", "encounter_rate", "count", "abundance")) %>%
    as.data.frame()
  
  # Store dataset sizes
  
  dataset_sizes$train[dataset_sizes$dataset == "all"] <- nrow(CBCH_all_train)
  dataset_sizes$test[dataset_sizes$dataset == "all"] <- nrow(CBCH_all_test)
  
  print(paste0("Bootstrap ", i, ": All Model Done (3/", 3+length(counties), ")"))
  
  ####################### MODEL 4: IN COUNTY MODELS ############################
  
  # We want to build a model based on data only from certain counties to see if
  # our pattern is replicable in subunits of our region.
  
  ### In Counties
  
  for(j in c("wa", "or", "mv")) {
    
    counties_select <- `if`(j == "wa", wa_counties, `if`(j == "or", or_counties, mv_counties))
    
    r_select <- `if`(j == "wa", r_wa, `if`(j == "or", r_or, r_mv))
    
    # Filter our training set to those points outside of the BCCH range.
    
    CBCH_county_train <- CBCH_ss[CBCH_ss$type == "train" & CBCH_ss$checklist_id %in% CBCH$checklist_id[CBCH$county %in% counties_select],] %>%
      select(names(CBCH_out_train))
    
    BCCH_train <- BCCH_ss[BCCH_ss$type == "train" & BCCH_ss$checklist_id %in% BCCH$checklist_id[BCCH$county %in% counties_select],] %>%
      select(names(BCCH_train))
    
    # First, we need to define a BCCH range.
    
    # Build step one of our two step hurdle model for BCCH - modeling encounter rate.
    
    # Calculate detection frequency for the balance random forest
    
    BCCH_detfreq <- mean(BCCH_train$species_observed)
    
    # Train a random forest model for encounter rate
    
    BCCH_train_er <- select(BCCH_train, -observation_count) %>%
      filter(is.na(time_observations_started) == FALSE)
    
    BCCH_er_model <- ranger(formula =  as.factor(species_observed) ~ ., 
                            data = BCCH_train_er,
                            importance = "permutation",
                            probability = TRUE,
                            replace = TRUE,
                            sample.fraction = c(BCCH_detfreq, BCCH_detfreq))
    
    # Select the mcc-f1 optimizing occurrence threshold
    
    BCCH_obs_pred <- tibble(obs = as.integer(BCCH_train_er$species_observed), 
                            pred = BCCH_er_model$predictions[, 2])
    
    BCCH_mcc_f1 <- mccf1(response = BCCH_obs_pred$obs, predictor = BCCH_obs_pred$pred)
    
    BCCH_threshold <- summary(BCCH_mcc_f1)$best_threshold[1]
    
    # Now, predict from our first model the BCCH abundance at each checklist in the
    # range of BCCH.
    
    # First, predict encounter rate
    
    BCCH_pred_er <- predict(BCCH_er_model,
                            CBCH_county_train,
                            type = "response")
    
    BCCH_pred_er <- as.data.frame(BCCH_pred_er$predictions[,2]) %>%
      setNames("predicted_er")
    
    # Apply threshold
    
    BCCH_pred_er$predicted_binary <- as.integer(BCCH_pred_er$predicted_er > BCCH_threshold)
    
    # Join to training data
    
    CBCH_county_train$BCCH_range_binary <- BCCH_pred_er$predicted_binary
    
    # Now that we have estimates of BCCH range at each point, build model (2)
    
    # Filter training data to cells with BCCH predicted to occur
    
    CBCH_county_polygon_train <- CBCH_county_train
    
    CBCH_county_train <- CBCH_county_train %>%
      filter(BCCH_range_binary == 1) %>%
      select(-BCCH_range, -BCCH_range_binary)
    
    # Build step one of our two step hurdle model - modeling encounter rate.
    
    # Calculate detection frequency for the balance random forest
    
    CBCH_county_detfreq <- mean(CBCH_county_train$species_observed)
    
    # Train a random forest model for encounter rate
    
    CBCH_county_train_er <- select(CBCH_county_train, -observation_count)
    
    CBCH_county_er_model <- ranger(formula =  as.factor(species_observed) ~ ., 
                                   data = CBCH_county_train_er,
                                   importance = "permutation",
                                   probability = TRUE,
                                   replace = TRUE,
                                   sample.fraction = c(CBCH_county_detfreq, CBCH_county_detfreq))
    
    # Select the mcc-f1 optimizing occurrence threshold
    
    CBCH_county_obs_pred <- tibble(obs = as.integer(CBCH_county_train_er$species_observed), 
                                   pred = CBCH_county_er_model$predictions[, 2])
    
    CBCH_county_mcc_f1 <- mccf1(response = CBCH_county_obs_pred$obs,
                                predictor = CBCH_county_obs_pred$pred)
    
    CBCH_county_threshold <- summary(CBCH_county_mcc_f1)$best_threshold[1]
    
    # Calibration model
    
    CBCH_county_calibration_model <- scam(obs ~ s(pred, k = 6, bs = "mpi"), 
                                          gamma = 2,
                                          data = CBCH_county_obs_pred)
    
    # Build step two of our two step hurdle model - modeling expected count.
    
    # Attach the predicted encounter rate based on in of bag samples
    
    CBCH_county_train_count <- CBCH_county_train
    
    CBCH_county_train_count$pred_er <- CBCH_county_er_model$predictions[, 2]
    
    # Subset to only observed or predicted detections
    
    CBCH_county_train_count <- CBCH_county_train_count %>% 
      filter(!is.na(observation_count),
             observation_count > 0 | pred_er > CBCH_county_threshold) %>% 
      select(-species_observed, -pred_er)
    
    # Model counts!
    
    CBCH_county_count_model <- ranger(formula = observation_count ~ .,
                                      data = CBCH_county_train_count,
                                      importance = "permutation",
                                      replace = TRUE,
                                      quantreg = TRUE, num.threads = cores)
    
    # Check model performance
    
    # Get the test set held in from training
    
    CBCH_county_test <- CBCH_ss[CBCH_ss$type == "test" & CBCH_ss$checklist_id %in% CBCH$checklist_id[CBCH$county %in% counties_select],] %>%
      select(names(CBCH_out_test))
    
    # Now, predict from our first model the BCCH binary occurence at each checklist 
    # in the range of BCCH.
    
    # First, predict encounter rate
    
    BCCH_county_pred_er <- predict(BCCH_er_model,
                                   CBCH_county_test,
                                   type = "response")
    
    BCCH_county_pred_er <- as.data.frame(BCCH_county_pred_er$predictions[,2]) %>%
      setNames("predicted_er")
    
    # Apply threshold
    
    BCCH_county_pred_er$predicted_binary <- as.integer(BCCH_county_pred_er$predicted_er > BCCH_threshold)
    
    # Join to training data
    
    CBCH_county_test$BCCH_range_binary <- BCCH_county_pred_er$predicted_binary
    
    # Filter training data to cells with BCCH predicted to occur
    
    CBCH_county_polygon_test <- CBCH_county_test
    
    CBCH_county_test <- CBCH_county_test %>%
      filter(BCCH_range_binary == 1) %>%
      select(-BCCH_range, -BCCH_range_binary)
    
    # Predict on test data using random forest model
    
    CBCH_county_testPred_er <- predict(CBCH_county_er_model, data = CBCH_county_test, type = "response")
    
    # Extract probability of detection
    
    CBCH_county_testPred_er <- CBCH_county_testPred_er$predictions[, 2]
    
    # Convert to binary using the threshold
    
    CBCH_county_testPred_binary <- as.integer(CBCH_county_testPred_er > CBCH_county_threshold)
    
    # Calibrate
    
    CBCH_county_testPred_calibrated <- predict(CBCH_county_calibration_model, 
                                               newdata = data.frame(pred = CBCH_county_testPred_er), 
                                               type = "response") %>% 
      as.numeric()
    
    # Add encounter rate to predict count
    
    CBCH_county_test$predicted_er <- CBCH_county_testPred_er
    
    # Estimate count
    
    CBCH_county_testPred_count <- predict(CBCH_county_count_model, data = CBCH_county_test, type = "response")
    
    CBCH_county_testPred_count <- CBCH_county_testPred_count$predictions
    
    # Relative abundance is the product of encounter rate and count
    
    CBCH_county_testPred_abundance <- CBCH_county_testPred_calibrated * CBCH_county_testPred_count
    
    # Combine all estimates together and join to actual observations
    
    CBCH_county_testPred_obs <- data.frame(id = seq_along(CBCH_county_testPred_abundance),
                                           # actual detection/non-detection
                                           obs_detection = as.integer(CBCH_county_test$species_observed),
                                           obs_count = CBCH_county_test$observation_count,
                                           # model estimates
                                           pred_count = CBCH_county_testPred_count,
                                           pred_abundance = CBCH_county_testPred_abundance,
                                           # binary detection/on-detection prediction
                                           pred_binary = CBCH_county_testPred_binary,
                                           # calibrated encounter rate
                                           pred_calibrated = CBCH_county_testPred_calibrated) %>%
      # constrain probabilities to 0-1
      mutate(pred_calibrated = pmin(pmax(pred_calibrated, 0), 1))
    
    # Mean squared error (MSE)
    
    CBCH_county_mse <- mean((CBCH_county_testPred_obs$obs_detection - CBCH_county_testPred_obs$pred_calibrated)^2, na.rm = TRUE)
    
    # spearman correlation, based on in range observations only
    
    CBCH_county_spearman_er <- cor(CBCH_county_testPred_obs$pred_calibrated[CBCH_county_testPred_obs$pred_binary > 0], 
                                   CBCH_county_testPred_obs$obs_detection[CBCH_county_testPred_obs$pred_binary > 0], 
                                   method = "spearman")
    
    # Precision-recall AUC
    
    CBCH_county_em <- precrec::evalmod(scores = CBCH_county_testPred_obs$pred_binary, 
                                       labels = CBCH_county_testPred_obs$obs_detection)
    
    CBCH_county_pr_auc <- precrec::auc(CBCH_county_em) %>% 
      filter(curvetypes == "PRC") %>% 
      pull(aucs)
    
    CBCH_county_roc_auc <- precrec::auc(CBCH_county_em) %>% 
      filter(curvetypes == "ROC") %>% 
      pull(aucs)
    
    # Calculate metrics for binary prediction: kappa, sensitivity, specificity
    
    CBCH_county_pa_metrics <- CBCH_county_testPred_obs %>% 
      select(id, obs_detection, pred_binary) %>% 
      presence.absence.accuracy(na.rm = TRUE, st.dev = FALSE)
    
    # MCC and F1
    
    CBCH_county_mcc_f1_test <- calculate_mcc_f1(CBCH_county_testPred_obs$obs_detection, 
                                                CBCH_county_testPred_obs$pred_binary)
    
    # subset to only those checklists where detection is predicted
    
    CBCH_county_testPred_detections <- filter(CBCH_county_testPred_obs, pred_binary > 0, 
                                              is.na(obs_count) == FALSE)
    
    # Count metrics
    
    CBCH_county_spearman_count <- cor(CBCH_county_testPred_detections$pred_count, 
                                      CBCH_county_testPred_detections$obs_count,
                                      method = "spearman")
    
    CBCH_county_pearson_log_count <- cor(log(CBCH_county_testPred_detections$pred_count + 1),
                                         log(CBCH_county_testPred_detections$obs_count + 1),
                                         method = "pearson")
    
    # Abundance metrics
    
    CBCH_county_spearman_abundance <- cor(CBCH_county_testPred_detections$pred_abundance, 
                                          CBCH_county_testPred_detections$obs_count,
                                          method = "spearman")
    
    CBCH_county_pearson_log_abundance <- cor(log(CBCH_county_testPred_detections$pred_abundance + 1),
                                             log(CBCH_county_testPred_detections$obs_count + 1),
                                             method = "pearson")
    
    # Combine metrics together
    
    counties_assessment[[j]][[i]] <- tibble(
      MSE = CBCH_county_mse,
      ER_Spearman = CBCH_county_spearman_er,
      Sensitivity = CBCH_county_pa_metrics$sensitivity,
      Specificity = CBCH_county_pa_metrics$specificity,
      Kappa = CBCH_county_pa_metrics$Kappa,
      PR_AUC = CBCH_county_pr_auc,
      ROC_AUC = CBCH_county_roc_auc,
      MCC = CBCH_county_mcc_f1_test$mcc,
      F1 = CBCH_county_mcc_f1_test$f1,
      Count_Spearman = CBCH_county_spearman_count,
      Log_Count_Pearson = CBCH_county_pearson_log_count,
      Abundance_Spearman = CBCH_county_spearman_abundance,
      Log_Abundance_Pearson = CBCH_county_pearson_log_abundance
    )
    
    # Now, predict!
    
    # Find peak time of day for CBCH detection from partial dependence
    
    CBCH_county_pd_time <- calculate_pd("time_observations_started",
                                        er_model = CBCH_county_er_model, 
                                        data = CBCH_county_train,
                                        calibration_model = CBCH_county_calibration_model,
                                        # make estimates at 30 minute intervals
                                        # using a subset of the training dataset
                                        x_res = 2 * 24, n = 1000) %>% 
      select(time_observations_started = x, encounter_rate)
    
    # hours with at least 1% of checklists
    
    CBCH_county_search_times <- CBCH_county_train %>% 
      mutate(time_observations_started = floor(time_observations_started)) %>%
      count(time_observations_started) %>% 
      mutate(pct = n / sum(n)) %>% 
      filter(pct >= 0.01)
    
    CBCH_county_t_peak <- CBCH_county_pd_time %>% 
      filter(floor(time_observations_started) %in% CBCH_county_search_times$time_observations_started) %>% 
      slice_max(order_by = encounter_rate) %>% 
      pull(time_observations_started)
    
    # Add effort covariates to prediction grid
    
    county_maximums <- pred_grid_counties %>%
      vect(geom = c("x", "y"), crs = crs, keepgeom = TRUE) %>%
      terra::intersect(as.polygons(r_select), .) %>%
      as.data.frame()
    
    # Filter prediction grid to values found within the modeled county.
    
    CBCH_county_pred_grid_eff <- pred_grid_counties %>%
      filter(road_dens <= max(county_maximums$road_dens), nighttime_light <= max(county_maximums$nighttime_light),
             pland_c13_urban <= max(county_maximums$pland_c13_urban), canopy_height <= max(county_maximums$canopy_height),
             pland_c01_evergreen_needleleaf <= max(county_maximums$pland_c01_evergreen_needleleaf)) %>%
      as.data.frame() %>%
      mutate(observation_date = ymd("2022-06-01"),
             year = year(observation_date),
             day_of_year = yday(observation_date),
             time_observations_started = CBCH_county_t_peak,
             duration_minutes = 120,
             effort_distance_km = 1,
             number_observers = 1)
    
    # Predict encounter rate first.
    
    CBCH_county_pred_er <- predict(CBCH_county_er_model, data = CBCH_county_pred_grid_eff, 
                                   type = "response")
    
    CBCH_county_pred_er <- as.data.frame(cbind(CBCH_county_pred_grid_eff$cell_id,
                                               CBCH_county_pred_er$predictions[,2])) %>%
      setNames(c("cell_id", "predicted_er"))
    
    # Apply threshold
    
    CBCH_county_pred_er$predicted_binary <- as.integer(CBCH_county_pred_er$predicted_er > CBCH_county_threshold)
    
    # Apply calibration
    
    CBCH_county_pred_er$calibrated_er <- predict(CBCH_county_calibration_model, 
                                                 data.frame(pred = CBCH_county_pred_er$predicted_er), 
                                                 type = "response") %>% 
      as.numeric()
    
    # Add predicted encounter rate required for count estimates
    
    CBCH_county_pred_grid_eff <- left_join(CBCH_county_pred_grid_eff, 
                                           CBCH_county_pred_er, 
                                           by = "cell_id")
    
    # Predict count estimate
    
    CBCH_county_pred_count <- predict(CBCH_county_count_model, 
                                      data = CBCH_county_pred_grid_eff, 
                                      type = "response")
    
    CBCH_county_pred_count <- as.data.frame(cbind(CBCH_county_pred_grid_eff$cell_id,
                                                  CBCH_county_pred_count$predictions)) %>%
      setNames(c("cell_id", "count"))
    
    
    # Add estimates to prediction grid
    
    CBCH_county_pred_grid_eff <- left_join(CBCH_county_pred_grid_eff, 
                                           CBCH_county_pred_count,
                                           by = "cell_id")
    
    # Bind predictions together
    
    CBCH_county_predictions <- CBCH_county_pred_grid_eff %>% 
      mutate(in_range = predicted_binary, encounter_rate = calibrated_er) %>%
      select(cell_id, x, y, in_range, encounter_rate, count) %>% 
      mutate(encounter_rate = pmin(pmax(encounter_rate, 0), 1))
    
    # Calculate relative abundance and rasterize.
    
    counties_r_pred[[j]][[i]] <- CBCH_county_predictions %>% 
      
      # Estimate relative abundance
      mutate(abundance = encounter_rate * count) %>% 
      
      # Convert to spatial features
      st_as_sf(coords = c("x", "y"), crs = crs) %>% 
      select(cell_id, in_range, encounter_rate, count, abundance) %>%
      
      # Rasterize
      rasterize(rast(r_select), field = c("cell_id", "in_range", "encounter_rate", "count", "abundance"),
                fun = "mean") %>% 
      setNames(c("cell_id", "in_range", "encounter_rate", "count", "abundance")) %>%
      as.data.frame()
  
    # Store dataset sizes
    
    dataset_sizes$train[dataset_sizes$dataset == j] <- nrow(CBCH_county_train)
    dataset_sizes$test[dataset_sizes$dataset == j] <- nrow(CBCH_county_test)
    
    print(paste0("Bootstrap ", i, ": ", str_to_upper(j) ," Model Done (", 3 + which(counties == j), "/", 3+length(counties), ")"))
    
  }
  
  ####################### MODEL 5: OUT COUNTY MODELS ############################
  
  for(j in c("vi", "sf", "nc")) {
    
    counties_select <- `if`(j == "vi", vi_counties, `if`(j == "sf", sf_counties, nc_counties))
    
    r_select <- `if`(j == "vi", r_vi, `if`(j == "sf", r_sf, r_nc))
    
    # Filter our training set to those points outside of the BCCH range.
    
    CBCH_county_train <- CBCH_ss[CBCH_ss$type == "train" & CBCH_ss$checklist_id %in% CBCH$checklist_id[CBCH$county %in% counties_select],] %>%
      select(names(CBCH_out_train))
    
    # Build step one of our two step hurdle model - modeling encounter rate.
    
    # Calculate detection frequency for the balance random forest
    
    CBCH_county_detfreq <- mean(CBCH_county_train$species_observed)
    
    # Train a random forest model for encounter rate
    
    CBCH_county_train_er <- select(CBCH_county_train, -observation_count, -BCCH_range)
    
    CBCH_county_er_model <- ranger(formula =  as.factor(species_observed) ~ ., 
                                   data = CBCH_county_train_er,
                                   importance = "permutation",
                                   probability = TRUE,
                                   replace = TRUE,
                                   sample.fraction = c(CBCH_county_detfreq, CBCH_county_detfreq))
    
    # Select the mcc-f1 optimizing occurrence threshold
    
    CBCH_county_obs_pred <- tibble(obs = as.integer(CBCH_county_train_er$species_observed), 
                                   pred = CBCH_county_er_model$predictions[, 2])
    
    CBCH_county_mcc_f1 <- mccf1(response = CBCH_county_obs_pred$obs,
                                predictor = CBCH_county_obs_pred$pred)
    
    CBCH_county_threshold <- summary(CBCH_county_mcc_f1)$best_threshold[1]
    
    # Calibration model
    
    CBCH_county_calibration_model <- scam(obs ~ s(pred, k = 6, bs = "mpi"), 
                                          gamma = 2,
                                          data = CBCH_county_obs_pred)
    
    # Build step two of our two step hurdle model - modeling expected count.
    
    # Attach the predicted encounter rate based on out of bag samples
    
    CBCH_county_train_count <- CBCH_county_train %>%
      select(-BCCH_range)
    
    CBCH_county_train_count$pred_er <- CBCH_county_er_model$predictions[, 2]
    
    # Subset to only observed or predicted detections
    
    CBCH_county_train_count <- CBCH_county_train_count %>% 
      filter(!is.na(observation_count),
             observation_count > 0 | pred_er > CBCH_county_threshold) %>% 
      select(-species_observed, -pred_er)
    
    # Model counts!
    
    CBCH_county_count_model <- ranger(formula = observation_count ~ .,
                                      data = CBCH_county_train_count,
                                      importance = "permutation",
                                      replace = TRUE,
                                      quantreg = TRUE, num.threads = cores)
    
    # Check model performance
    
    # Get the test set held out from training
    
    CBCH_county_test <- CBCH_ss[CBCH_ss$type == "test" & CBCH_ss$checklist_id %in% CBCH$checklist_id[CBCH$county %in% counties_select],] %>%
      select(names(CBCH_out_train))
    
    # Predict on test data using random forest model
    
    CBCH_county_testPred_er <- predict(CBCH_county_er_model, data = CBCH_county_test, type = "response")
    
    # Extract probability of detection
    
    CBCH_county_testPred_er <- CBCH_county_testPred_er$predictions[, 2]
    
    # Convert to binary using the threshold
    
    CBCH_county_testPred_binary <- as.integer(CBCH_county_testPred_er > CBCH_county_threshold)
    
    # Calibrate
    
    CBCH_county_testPred_calibrated <- predict(CBCH_county_calibration_model, 
                                               newdata = data.frame(pred = CBCH_county_testPred_er), 
                                               type = "response") %>% 
      as.numeric()
    
    # Add encounter rate to predict count
    
    CBCH_county_test$predicted_er <- CBCH_county_testPred_er
    
    # Estimate count
    
    CBCH_county_testPred_count <- predict(CBCH_county_count_model, data = CBCH_county_test, type = "response")
    
    CBCH_county_testPred_count <- CBCH_county_testPred_count$predictions
    
    # Relative abundance is the product of encounter rate and count
    
    CBCH_county_testPred_abundance <- CBCH_county_testPred_calibrated * CBCH_county_testPred_count
    
    # Combine all estimates together and join to actual observations
    
    CBCH_county_testPred_obs <- data.frame(id = seq_along(CBCH_county_testPred_abundance),
                                           # actual detection/non-detection
                                           obs_detection = as.integer(CBCH_county_test$species_observed),
                                           obs_count = CBCH_county_test$observation_count,
                                           # model estimates
                                           pred_count = CBCH_county_testPred_count,
                                           pred_abundance = CBCH_county_testPred_abundance,
                                           # binary detection/on-detection prediction
                                           pred_binary = CBCH_county_testPred_binary,
                                           # calibrated encounter rate
                                           pred_calibrated = CBCH_county_testPred_calibrated) %>%
      # constrain probabilities to 0-1
      mutate(pred_calibrated = pmin(pmax(pred_calibrated, 0), 1))
    
    # Mean squared error (MSE)
    
    CBCH_county_mse <- mean((CBCH_county_testPred_obs$obs_detection - CBCH_county_testPred_obs$pred_calibrated)^2, na.rm = TRUE)
    
    # spearman correlation, based on in range observations only
    
    CBCH_county_spearman_er <- cor(CBCH_county_testPred_obs$pred_calibrated[CBCH_county_testPred_obs$pred_binary > 0], 
                                   CBCH_county_testPred_obs$obs_detection[CBCH_county_testPred_obs$pred_binary > 0], 
                                   method = "spearman")
    
    # Precision-recall AUC
    
    CBCH_county_em <- precrec::evalmod(scores = CBCH_county_testPred_obs$pred_binary, 
                                       labels = CBCH_county_testPred_obs$obs_detection)
    
    CBCH_county_pr_auc <- precrec::auc(CBCH_county_em) %>% 
      filter(curvetypes == "PRC") %>% 
      pull(aucs)
    
    CBCH_county_roc_auc <- precrec::auc(CBCH_county_em) %>% 
      filter(curvetypes == "ROC") %>% 
      pull(aucs)
    
    # Calculate metrics for binary prediction: kappa, sensitivity, specificity
    
    CBCH_county_pa_metrics <- CBCH_county_testPred_obs %>% 
      select(id, obs_detection, pred_binary) %>% 
      presence.absence.accuracy(na.rm = TRUE, st.dev = FALSE)
    
    # MCC and F1
    
    CBCH_county_mcc_f1_test <- calculate_mcc_f1(CBCH_county_testPred_obs$obs_detection, 
                                                CBCH_county_testPred_obs$pred_binary)
    
    # subset to only those checklists where detection is predicted
    
    CBCH_county_testPred_detections <- filter(CBCH_county_testPred_obs, pred_binary > 0, 
                                              is.na(obs_count) == FALSE)
    
    # Count metrics
    
    CBCH_county_spearman_count <- cor(CBCH_county_testPred_detections$pred_count, 
                                      CBCH_county_testPred_detections$obs_count,
                                      method = "spearman")
    
    CBCH_county_pearson_log_count <- cor(log(CBCH_county_testPred_detections$pred_count + 1),
                                         log(CBCH_county_testPred_detections$obs_count + 1),
                                         method = "pearson")
    
    # Abundance metrics
    
    CBCH_county_spearman_abundance <- cor(CBCH_county_testPred_detections$pred_abundance, 
                                          CBCH_county_testPred_detections$obs_count,
                                          method = "spearman")
    
    CBCH_county_pearson_log_abundance <- cor(log(CBCH_county_testPred_detections$pred_abundance + 1),
                                             log(CBCH_county_testPred_detections$obs_count + 1),
                                             method = "pearson")
    
    # Combine metrics together
    
    counties_assessment[[j]][[i]] <- tibble(
      MSE = CBCH_county_mse,
      ER_Spearman = CBCH_county_spearman_er,
      Sensitivity = CBCH_county_pa_metrics$sensitivity,
      Specificity = CBCH_county_pa_metrics$specificity,
      Kappa = CBCH_county_pa_metrics$Kappa,
      PR_AUC = CBCH_county_pr_auc,
      ROC_AUC = CBCH_county_roc_auc,
      MCC = CBCH_county_mcc_f1_test$mcc,
      F1 = CBCH_county_mcc_f1_test$f1,
      Count_Spearman = CBCH_county_spearman_count,
      Log_Count_Pearson = CBCH_county_pearson_log_count,
      Abundance_Spearman = CBCH_county_spearman_abundance,
      Log_Abundance_Pearson = CBCH_county_pearson_log_abundance
    )
    
    # Now, predict!
    
    # Find peak time of day for CBCH detection from partial dependence
    
    CBCH_county_pd_time <- calculate_pd("time_observations_started",
                                        er_model = CBCH_county_er_model, 
                                        data = CBCH_county_train,
                                        calibration_model = CBCH_county_calibration_model,
                                        # make estimates at 30 minute intervals
                                        # using a subset of the training dataset
                                        x_res = 2 * 24, n = 1000) %>% 
      select(time_observations_started = x, encounter_rate)
    
    # hours with at least 1% of checklists
    
    CBCH_county_search_times <- CBCH_county_train %>% 
      mutate(time_observations_started = floor(time_observations_started)) %>%
      count(time_observations_started) %>% 
      mutate(pct = n / sum(n)) %>% 
      filter(pct >= 0.01)
    
    CBCH_county_t_peak <- CBCH_county_pd_time %>% 
      filter(floor(time_observations_started) %in% CBCH_county_search_times$time_observations_started) %>% 
      slice_max(order_by = encounter_rate) %>% 
      pull(time_observations_started)
    
    # Add effort covariates to prediction grid
  
    county_maximums <- pred_grid_counties %>%
      vect(geom = c("x", "y"), crs = crs, keepgeom = TRUE) %>%
      terra::intersect(as.polygons(r_select), .) %>%
      as.data.frame()
    
    # Filter prediction grid to values found within the modeled county.
    
    CBCH_county_pred_grid_eff <- pred_grid_counties %>%
      filter(road_dens <= max(county_maximums$road_dens), nighttime_light <= max(county_maximums$nighttime_light),
             pland_c13_urban <= max(county_maximums$pland_c13_urban), canopy_height <= max(county_maximums$canopy_height),
             pland_c01_evergreen_needleleaf <= max(county_maximums$pland_c01_evergreen_needleleaf)) %>%
      as.data.frame() %>%
      mutate(observation_date = ymd("2022-06-01"),
             year = year(observation_date),
             day_of_year = yday(observation_date),
             time_observations_started = CBCH_county_t_peak,
             duration_minutes = 120,
             effort_distance_km = 1,
             number_observers = 1)
    
    # Predict encounter rate first.
    
    CBCH_county_pred_er <- predict(CBCH_county_er_model, data = CBCH_county_pred_grid_eff, 
                                   type = "response")
    
    CBCH_county_pred_er <- as.data.frame(cbind(CBCH_county_pred_grid_eff$cell_id,
                                               CBCH_county_pred_er$predictions[,2])) %>%
      setNames(c("cell_id", "predicted_er"))
    
    # Apply threshold
    
    CBCH_county_pred_er$predicted_binary <- as.integer(CBCH_county_pred_er$predicted_er > CBCH_county_threshold)
    
    # Apply calibration
    
    CBCH_county_pred_er$calibrated_er <- predict(CBCH_county_calibration_model, 
                                                 data.frame(pred = CBCH_county_pred_er$predicted_er), 
                                                 type = "response") %>% 
      as.numeric()
    
    # Add predicted encounter rate required for count estimates
    
    CBCH_county_pred_grid_eff <- left_join(CBCH_county_pred_grid_eff, 
                                           CBCH_county_pred_er, 
                                           by = "cell_id")
    
    # Predict count estimate
    
    CBCH_county_pred_count <- predict(CBCH_county_count_model, 
                                      data = CBCH_county_pred_grid_eff, 
                                      type = "response")
    
    CBCH_county_pred_count <- as.data.frame(cbind(CBCH_county_pred_grid_eff$cell_id,
                                                  CBCH_county_pred_count$predictions)) %>%
      setNames(c("cell_id", "count"))
    
    
    # Add estimates to prediction grid
    
    CBCH_county_pred_grid_eff <- left_join(CBCH_county_pred_grid_eff, 
                                           CBCH_county_pred_count,
                                           by = "cell_id")
    
    # Bind predictions together
    
    CBCH_county_predictions <- CBCH_county_pred_grid_eff %>% 
      mutate(in_range = predicted_binary, encounter_rate = calibrated_er) %>%
      select(cell_id, x, y, in_range, encounter_rate, count) %>% 
      mutate(encounter_rate = pmin(pmax(encounter_rate, 0), 1))
    
    # Calculate relative abundance and rasterize.
    
    counties_r_pred[[j]][[i]] <- CBCH_county_predictions %>% 
      
      # Estimate relative abundance
      mutate(abundance = encounter_rate * count) %>% 
      
      # Convert to spatial features
      st_as_sf(coords = c("x", "y"), crs = crs) %>% 
      select(cell_id, in_range, encounter_rate, count, abundance) %>%
      
      # Rasterize
      rasterize(rast(r_select), field = c("cell_id", "in_range", "encounter_rate", "count", "abundance"),
                fun = "mean") %>% 
      setNames(c("cell_id", "in_range", "encounter_rate", "count", "abundance")) %>%
      as.data.frame()
    
    # Store dataset sizes
    
    dataset_sizes$train[dataset_sizes$dataset == j] <- nrow(CBCH_county_train)
    dataset_sizes$test[dataset_sizes$dataset == j] <- nrow(CBCH_county_test)
    
    print(paste0("Bootstrap ", i, ": ", str_to_upper(j) ," Model Done (", 3 + which(counties == j), "/", 3+length(counties), ")"))
    
  }
  
  ############################## DATASET SIZES #################################
  
  n[[i]] <- dataset_sizes
  
  
  # Periodic save of results
  
  if(i %in% seq(from = 0, to = 100, by = 5)) {
    
    print(paste0("Saving ", i, " bootstraps"))
    
    save.image(save_path)
    
  }
  
  
}

  ############################ POST-LOOP ANALYSES ###############################



############################# SAMPLE SIZES #####################################

n_summary <- data.frame(loc = rep(c("out", "all", "in", counties), each = 2),
                        type = rep(c("train", "test"), times = 3 + length(counties)))

for(i in seq_along(seeds)[1:length(n)]) {
  
  for(j in n_summary$loc) {
    
    for(k in n_summary$type) {
      
      n_summary[n_summary$loc == j & n_summary$type == k, paste0("boot", seeds[i])] <- ifelse(k == "train",
                                                                                              n[[i]]$train[n[[i]]$dataset == j],
                                                                                              n[[i]]$test[n[[i]]$dataset == j])
      
    }
  }
  
}

n_summary$n_mean <- NA
n_summary$n_sd <- NA

for(i in 1:nrow(n_summary)) {
  
  n_summary$n_mean[i] <- mean(unname(unlist(n_summary[i, paste0("boot", seeds[1:length(n)])])))
  n_summary$n_sd[i] <- sd(unname(unlist(n_summary[i, paste0("boot", seeds[1:length(n)])])))
  
}


write_csv(n_summary[,c("loc", "type", "n_mean", "n_sd")], "./Outputs/Tables/n_table.csv")

n_table <- kable(n_summary[,c("loc", "type", "n_mean", "n_sd")])


############################ PLOTTING PREDICTIONS ##############################

abundance_predictions <- data.frame(cell_id = rep(out_r_pred[[1]]$cell_id, times = length(locs[1:1])),
                                    loc = rep(locs[1:3], each = length(out_r_pred[[1]]$cell_id)))

abundance_predictions_counties <- data.frame(cell_id = c(counties_r_pred[["wa"]][[1]]$cell_id, counties_r_pred[["or"]][[1]]$cell_id,
                                                         counties_r_pred[["mv"]][[1]]$cell_id, counties_r_pred[["vi"]][[1]]$cell_id,
                                                         counties_r_pred[["sf"]][[1]]$cell_id, counties_r_pred[["nc"]][[1]]$cell_id),
                                             loc = rep(locs[4:length(locs)], times = c(length(counties_r_pred[["wa"]][[1]]$cell_id), length(counties_r_pred[["or"]][[1]]$cell_id),
                                                                                      length(counties_r_pred[["mv"]][[1]]$cell_id), length(counties_r_pred[["vi"]][[1]]$cell_id),
                                                                                      length(counties_r_pred[["sf"]][[1]]$cell_id), length(counties_r_pred[["nc"]][[1]]$cell_id))))

for(j in seq_along(seeds)) {
  
  out_predictions <- out_r_pred[[j]]$abundance * out_r_pred[[j]]$in_range
  abundance_predictions[abundance_predictions$loc == "Out", paste0("boot", seeds[j])] <- out_predictions
  
  in_predictions <- in_r_pred[[j]]$abundance * in_r_pred[[j]]$in_range
  abundance_predictions[abundance_predictions$loc == "In", paste0("boot", seeds[j])] <- in_predictions
  
  all_predictions <- all_r_pred[[j]]$abundance * all_r_pred[[j]]$in_range
  abundance_predictions[abundance_predictions$loc == "All", paste0("boot", seeds[j])] <- all_predictions
  
  wa_predictions <- counties_r_pred[["wa"]][[j]]$abundance * counties_r_pred[["wa"]][[j]]$in_range
  abundance_predictions_counties[abundance_predictions_counties$loc == "WA", paste0("boot", seeds[j])] <- wa_predictions
  
  or_predictions <- counties_r_pred[["or"]][[j]]$abundance * counties_r_pred[["or"]][[j]]$in_range
  abundance_predictions_counties[abundance_predictions_counties$loc == "OR", paste0("boot", seeds[j])] <- or_predictions
  
  mv_predictions <- counties_r_pred[["mv"]][[j]]$abundance * counties_r_pred[["mv"]][[j]]$in_range
  abundance_predictions_counties[abundance_predictions_counties$loc == "MV", paste0("boot", seeds[j])] <- mv_predictions
  
  vi_predictions <- counties_r_pred[["vi"]][[j]]$abundance * counties_r_pred[["vi"]][[j]]$in_range
  abundance_predictions_counties[abundance_predictions_counties$loc == "VI", paste0("boot", seeds[j])] <- vi_predictions

  sf_predictions <- counties_r_pred[["sf"]][[j]]$abundance * counties_r_pred[["sf"]][[j]]$in_range
  abundance_predictions_counties[abundance_predictions_counties$loc == "SF", paste0("boot", seeds[j])] <- sf_predictions
  
  nc_predictions <- counties_r_pred[["nc"]][[j]]$abundance * counties_r_pred[["nc"]][[j]]$in_range
  abundance_predictions_counties[abundance_predictions_counties$loc == "NC", paste0("boot", seeds[j])] <- nc_predictions
  
}

for(i in 1:nrow(abundance_predictions)) {
  
  abundance_predictions$median[i] <- median(unlist(c(unname(abundance_predictions[i,paste0("boot", seeds)]))))
  abundance_predictions$lwr[i] <- quantile(unlist(c(unname(abundance_predictions[i,paste0("boot", seeds)]))), 0.025)
  abundance_predictions$upr[i] <- quantile(unlist(c(unname(abundance_predictions[i,paste0("boot", seeds)]))), 0.975)
  
}

for(i in 1:nrow(abundance_predictions_counties)) {
  
  abundance_predictions_counties$median[i] <- median(unlist(c(unname(abundance_predictions_counties[i,paste0("boot", seeds)]))))
  abundance_predictions_counties$lwr[i] <- quantile(unlist(c(unname(abundance_predictions_counties[i,paste0("boot", seeds)]))), 0.025)
  abundance_predictions_counties$upr[i] <- quantile(unlist(c(unname(abundance_predictions_counties[i,paste0("boot", seeds)]))), 0.975)
  
}

abundance_predictions <- left_join(abundance_predictions, pred_grid, by = "cell_id") %>%
  mutate(loc = case_when(loc == "Out" ~ "Allopatric",
                         loc == "All" ~ "Region-wide",
                         loc == "In" ~ "Sympatric"))

abundance_predictions_counties <- left_join(abundance_predictions_counties, pred_grid, by = "cell_id") %>%
  mutate(loc_InOut = case_when(loc %in% c("VI", "SF", "NC") ~ "Allopatric",
                               loc %in% c("WA", "OR", "MV") ~ "Sympatric"))

# Nighttime Light

abundance_predictions <- abundance_predictions[abundance_predictions$nighttime_light <= 65,]
abundance_predictions_counties <- abundance_predictions_counties[abundance_predictions_counties$nighttime_light <= 65,]

# Road Density

abundance_predictions <- abundance_predictions[abundance_predictions$road_dens <= quantile(CBCH$road_dens[CBCH$BCCH_range == "Out"], 0.999),]
abundance_predictions_counties <- abundance_predictions_counties[abundance_predictions_counties$road_dens <= quantile(CBCH$road_dens[CBCH$BCCH_range == "Out"], 0.999),]

# Convert text columns to factors in order desired in plots.

abundance_predictions$loc <- factor(abundance_predictions$loc, levels = c("Allopatric", "Region-wide", "Sympatric"))
abundance_predictions_counties$loc_InOut <- factor(abundance_predictions_counties$loc_InOut, levels = c("Allopatric", "Sympatric"))

# Convert to long-form so we can fit a GAM smoother to each bootstrap

abund_long <- pivot_longer(abundance_predictions, cols = starts_with("boot"), names_to = "bootstrap", values_to = "prediction")
abund_long$boot_loc <- paste0(abund_long$bootstrap, abund_long$loc)

# Add a column identifying whether a cell is in a Sympatric (In) or Allopatric (Out) area.

abund_long <- mutate(abund_long, region = case_when(cell_id %in% in_grid$cell_id ~ "In",
                                        cell_id %in% out_grid$cell_id ~ "Out"))

# Subset to create sets we'll use for plots - "internal" sets refer to predictions
# from when a model is predicted back onto the area the data it was built from
# came from, and "external" sets refer to predictions from when a model is
# predicted outside of the area the data it was build from came from.

# Create 'internal' for use in Figure 2 - demonstrating the relationship to focal
# variables within the area a model was built on.

internal <- filter(abund_long, loc %in% as.factor("Allopatric") & region == "Out" | loc %in% as.factor("Sympatric") & region == "In")

# Create 'external_out' for use in Figures 5 (urban) and 6 (forest) - comparing the real-world predicted
# abundance and the predicted abundance under imagined scenarios in allopatric areas 

external_out <- abund_long %>%
  filter(loc %in% as.factor("Allopatric") & region == "Out" | loc %in% as.factor("Sympatric") & region == "Out") %>%
  mutate(scenario = case_when(loc == "Allopatric" ~ "Real-world",
                            loc == "Sympatric" ~ "BCCH Present"))

external_out$scenario = factor(external_out$scenario, c("Real-world", "BCCH Present"))

# Create 'external_in' for use in Figure 5 - comparing the real-world predicted
# abundance and the predicted abundance under imagined scenarios in sympatric areas

external_in <- abund_long %>%
  filter(loc %in% as.factor("Sympatric") & region == "In" | loc %in% as.factor("Allopatric") & region == "In") %>%
  mutate(scenario = case_when(loc == "Sympatric" ~ "Real-world",
                              loc == "Allopatric" ~ "BCCH Absent"))

external_in$scenario = factor(external_in$scenario, c("Real-world", "BCCH Absent"))

# Create plots for Road Density

abundPlot_Internal_RoadDensity <- ggplot() +
  scale_color_manual(values = plasma(3, begin = 0, end = 0.75, direction = -1)) +
  labs(x = "\n Road Density (m/km^2)",  y = "Predicted CBCH Rel. Abund. \n", colour = "Region") +
  geom_point(data = internal, aes(x = road_dens, y = prediction, colour = loc), alpha = 0.01, size = 2, show.legend = FALSE) +
  geom_smooth(data = internal, aes(x = road_dens, y = prediction, group = boot_loc), colour = "grey", alpha = 0.5, method = "gam", se = T) +
  geom_smooth(data = internal, aes(x = road_dens, y = prediction, colour = loc), method = "gam", se = T, size = 3) +
  scale_x_continuous(expand = c(0, 0), breaks = c(5000, 10000, 15000)) + 
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +
  theme(axis.text = element_text(size = 20), axis.title.x = element_text(size = 24, vjust = 2),
        axis.title.y = element_text(size = 24),
        legend.key.size = unit(1, "in"), legend.text = element_text(size = 20),
        legend.title = element_text(size = 20), axis.ticks.length = unit(0.1, "in")) +
  coord_cartesian(ylim = c(0, 4.2), xlim = c(0,16000))

abundPlot_ExternalOut_RoadDensity <- ggplot() +
  scale_color_manual(values = plasma(3, begin = 0, end = 0.75, direction = -1)) +
  labs(x = "\n Road Density (m/km^2)", y = "Predicted CBCH Rel. Abund. \n", colour = "Allopatric Areas") +
  geom_point(data = external_out, aes(x = road_dens, y = prediction, colour = scenario), alpha = 0.01, size = 2, show.legend = FALSE) +
  geom_smooth(data = external_out, aes(x = road_dens, y = prediction, group = boot_loc), colour = "grey", alpha = 0.5, method = "gam", se = T) +
  geom_smooth(data = external_out, aes(x = road_dens, y = prediction, colour = scenario), method = "gam", se = T, size = 3) +
  scale_x_continuous(expand = c(0, 0), breaks = c(5000, 10000, 15000)) + 
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +
  theme(axis.text = element_text(size = 20), axis.title.x = element_text(size = 16, vjust = 2),
        axis.title.y = element_text(size = 16),
        legend.key.size = unit(1, "in"), legend.text = element_text(size = 18),
        legend.title = element_text(size = 20, hjust = 0.5), axis.ticks.length = unit(0.1, "in"),
        legend.position = "top", legend.title.position = "top") +
  coord_cartesian(ylim = c(0, 4.2), xlim = c(0,16000))

abundPlot_ExternalIn_RoadDensity <- ggplot() +
  scale_color_manual(values = plasma(3, begin = 0, end = 0.75, direction = -1)) +
  labs(x = "\n Road Density (m/km^2)", colour = "Sympatric Areas") +
  geom_point(data = external_in, aes(x = road_dens, y = prediction, colour = scenario), alpha = 0.01, size = 2, show.legend = FALSE) +
  geom_smooth(data = external_in, aes(x = road_dens, y = prediction, group = boot_loc), colour = "grey", alpha = 0.5, method = "gam", se = T) +
  geom_smooth(data = external_in, aes(x = road_dens, y = prediction, colour = scenario), method = "gam", se = T, size = 3) +
  scale_x_continuous(expand = c(0, 0), breaks = c(5000, 10000, 15000)) + 
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +
  theme(axis.text = element_text(size = 20), axis.title.x = element_text(size = 16, vjust = 2),
        axis.title.y = element_blank(), axis.text.y = element_blank(),
        legend.key.size = unit(1, "in"), legend.text = element_text(size = 18),
        legend.title = element_text(size = 20, hjust = 0.5), axis.ticks.length = unit(0.1, "in"),
        legend.position = "top", legend.title.position = "top") +
  coord_cartesian(ylim = c(0, 4.2), xlim = c(0,16000))

# Create plots for Nighttime Light

abundPlot_Internal_NighttimeLight <- ggplot() +
  scale_color_manual(values = plasma(3, begin = 0, end = 0.75, direction = -1)) +
  labs(x = "\n Nighttime Light Radiance", colour = "Region") +
  geom_point(data = internal, aes(x = nighttime_light, y = prediction, colour = loc), alpha = 0.01, size = 2, show.legend = FALSE) +
  geom_smooth(data = internal, aes(x = nighttime_light, y = prediction, group = boot_loc), colour = "grey", alpha = 0.5, method = "gam", se = T) +
  geom_smooth(data = internal, aes(x = nighttime_light, y = prediction, colour = loc), method = "gam", se = T, size = 3) +
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +
  theme(axis.text = element_text(size = 20), axis.title.x = element_text(size = 24, vjust = 2),
        axis.text.y = element_blank(), axis.title.y = element_blank(),
        legend.key.size = unit(1, "in"), legend.text = element_text(size = 20),
        legend.title = element_text(size = 20), axis.ticks.length = unit(0.1, "in")) +
  coord_cartesian(ylim = c(0, 4.2), xlim = c(0,35))

abundPlot_ExternalOut_NighttimeLight <- ggplot() +
  scale_color_manual(values = plasma(3, begin = 0, end = 0.75, direction = -1)) +
  labs(x = "\n Nighttime Light Radiance", y = "Predicted CBCH Rel. Abund. \n", colour = "Scenario") +
  geom_point(data = external_out, aes(x = nighttime_light, y = prediction, colour = scenario), alpha = 0.01, size = 2, show.legend = FALSE) +
  geom_smooth(data = external_out, aes(x = nighttime_light, y = prediction, group = boot_loc), colour = "grey", alpha = 0.5, method = "gam", se = T) +
  geom_smooth(data = external_out, aes(x = nighttime_light, y = prediction, colour = scenario), method = "gam", se = T, size = 3) +
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +
  theme(axis.text = element_text(size = 20), axis.title.x = element_text(size = 16, vjust = 2),
        axis.title.y = element_text(size = 16),
        legend.key.size = unit(1, "in"), legend.text = element_text(size = 20),
        legend.title = element_text(size = 20, hjust = 0.5), axis.ticks.length = unit(0.1, "in"),
        legend.position = "none") +
  coord_cartesian(ylim = c(0, 4.2), xlim = c(0,35))

abundPlot_ExternalIn_NighttimeLight <- ggplot() +
  scale_color_manual(values = plasma(3, begin = 0, end = 0.75, direction = -1)) +
  labs(x = "\n Nighttime Light Radiance", colour = "Scenario") +
  geom_point(data = external_in, aes(x = nighttime_light, y = prediction, colour = scenario), alpha = 0.01, size = 2, show.legend = FALSE) +
  geom_smooth(data = external_in, aes(x = nighttime_light, y = prediction, group = boot_loc), colour = "grey", alpha = 0.5, method = "gam", se = T) +
  geom_smooth(data = external_in, aes(x = nighttime_light, y = prediction, colour = scenario), method = "gam", se = T, size = 3) +
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +
  theme(axis.text = element_text(size = 20), axis.title.x = element_text(size = 16, vjust = 2),
        axis.text.y = element_blank(), axis.title.y = element_blank(),
        legend.key.size = unit(1, "in"), legend.text = element_text(size = 20),
        legend.title = element_text(size = 20, hjust = 0.5), axis.ticks.length = unit(0.1, "in"),
        legend.position = "none") +
  coord_cartesian(ylim = c(0, 4.2), xlim = c(0,35))

# Create plots for Urban Cover

abundPlot_Internal_UrbanCover <- ggplot() +
  scale_color_manual(values = plasma(3, begin = 0, end = 0.75, direction = -1)) +
  labs(x = "\n Urban Land Cover (%)", colour = "Region") +
  geom_point(data = internal, aes(x = pland_c13_urban, y = prediction, colour = loc), alpha = 0.01, size = 2, show.legend = FALSE) +
  geom_smooth(data = internal, aes(x = pland_c13_urban, y = prediction, group = boot_loc), colour = "grey", alpha = 0.5, method = "gam", se = T) +
  geom_smooth(data = internal, aes(x = pland_c13_urban, y = prediction, colour = loc), method = "gam", se = T, size = 3) +
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +
  theme(axis.text = element_text(size = 20), axis.title.x = element_text(size = 24, vjust = 2),
        axis.text.y = element_blank(), axis.title.y = element_blank(),
        legend.key.size = unit(1, "in"), legend.text = element_text(size = 20),
        legend.title = element_text(size = 20), axis.ticks.length = unit(0.1, "in")) +
  coord_cartesian(ylim = c(0, 4.2), xlim = c(0,105))

abundPlot_ExternalOut_UrbanCover <- ggplot() +
  scale_color_manual(values = plasma(3, begin = 0, end = 0.75, direction = -1)) +
  labs(x = "\n Urban Land Cover (%)", y = "Predicted CBCH Rel. Abund. \n", colour = "Scenario") +
  geom_point(data = external_out, aes(x = pland_c13_urban, y = prediction, colour = scenario), alpha = 0.01, size = 2, show.legend = FALSE) +
  geom_smooth(data = external_out, aes(x = pland_c13_urban, y = prediction, group = boot_loc), colour = "grey", alpha = 0.5, method = "gam", se = T) +
  geom_smooth(data = external_out, aes(x = pland_c13_urban, y = prediction, colour = scenario), method = "gam", se = T, size = 3) +
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +
  theme(axis.text = element_text(size = 20), axis.title.x = element_text(size = 16, vjust = 2),
        axis.title.y = element_text(size = 16),
        legend.key.size = unit(1, "in"), legend.text = element_text(size = 20),
        legend.title = element_text(size = 20), axis.ticks.length = unit(0.1, "in"),
        legend.position = "none") +
  coord_cartesian(ylim = c(0, 4.2), xlim = c(0,105))

abundPlot_ExternalIn_UrbanCover <- ggplot() +
  scale_color_manual(values = plasma(3, begin = 0, end = 0.75, direction = -1)) +
  labs(x = "\n Urban Land Cover (%)", colour = "Scenario") +
  geom_point(data = external_in, aes(x = pland_c13_urban, y = prediction, colour = scenario), alpha = 0.01, size = 2, show.legend = FALSE) +
  geom_smooth(data = external_in, aes(x = pland_c13_urban, y = prediction, group = boot_loc), colour = "grey", alpha = 0.5, method = "gam", se = T) +
  geom_smooth(data = external_in, aes(x = pland_c13_urban, y = prediction, colour = scenario), method = "gam", se = T, size = 3) +
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +
  theme(axis.text = element_text(size = 20), axis.title.x = element_text(size = 16, vjust = 2),
        axis.text.y = element_blank(), axis.title.y = element_blank(),
        legend.key.size = unit(1, "in"), legend.text = element_text(size = 20),
        legend.title = element_text(size = 20), axis.ticks.length = unit(0.1, "in"),
        legend.position = "none") +
  coord_cartesian(ylim = c(0, 4.2), xlim = c(0,105))

# Create plots for Canopy Height

abundPlot_Internal_CanopyHeight <- ggplot() +
  scale_color_manual(values = plasma(3, begin = 0, end = 0.75, direction = -1)) +
  labs(x = "\n Mean Canopy Height (m)",  y = "Predicted CBCH Rel. Abund. \n", colour = "Region") +
  geom_point(data = internal, aes(x = canopy_height, y = prediction, colour = loc), alpha = 0.01, size = 2, show.legend = FALSE) +
  geom_smooth(data = internal, aes(x = canopy_height, y = prediction, group = boot_loc), colour = "grey", alpha = 0.5, method = "gam", se = T) +
  geom_smooth(data = internal, aes(x = canopy_height, y = prediction, colour = loc), method = "gam", se = T, size = 3) +
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +
  theme(axis.text = element_text(size = 20), axis.title.x = element_text(size = 24, vjust = 2),
        axis.title.y = element_text(size = 24),
        legend.key.size = unit(1, "in"), legend.text = element_text(size = 20),
        legend.title = element_text(size = 20), axis.ticks.length = unit(0.1, "in")) +
  coord_cartesian(ylim = c(0, 4.2), xlim = c(0,21))

abundPlot_ExternalOut_CanopyHeight <- ggplot() +
  scale_color_manual(values = plasma(3, begin = 0, end = 0.75, direction = -1)) +
  labs(x = "\n Mean Canopy Height (m)", y = "Predicted CBCH Rel. Abund. \n", colour = "Allopatric Areas") +
  geom_point(data = external_out, aes(x = canopy_height, y = prediction, colour = scenario), alpha = 0.01, size = 2, show.legend = FALSE) +
  geom_smooth(data = external_out, aes(x = canopy_height, y = prediction, group = boot_loc), colour = "grey", alpha = 0.5, method = "gam", se = T) +
  geom_smooth(data = external_out, aes(x = canopy_height, y = prediction, colour = scenario), method = "gam", se = T, size = 3) +
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +
  theme(axis.text = element_text(size = 20), axis.title.x = element_text(size = 16, vjust = 2),
        axis.title.y = element_text(size = 16),
        legend.key.size = unit(1, "in"), legend.text = element_text(size = 18),
        legend.title = element_text(size = 20, hjust = 0.5), axis.ticks.length = unit(0.1, "in"),
        legend.position = "top", legend.title.position = "top") +
  coord_cartesian(ylim = c(0, 4.2), xlim = c(0,21))

abundPlot_ExternalIn_CanopyHeight <- ggplot() +
  scale_color_manual(values = plasma(3, begin = 0, end = 0.75, direction = -1)) +
  labs(x = "\n Mean Canopy Height (m)", colour = "Sympatric Areas") +
  geom_point(data = external_in, aes(x = canopy_height, y = prediction, colour = scenario), alpha = 0.01, size = 2, show.legend = FALSE) +
  geom_smooth(data = external_in, aes(x = canopy_height, y = prediction, group = boot_loc), colour = "grey", alpha = 0.5, method = "gam", se = T) +
  geom_smooth(data = external_in, aes(x = canopy_height, y = prediction, colour = scenario), method = "gam", se = T, size = 3) +
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +
  theme(axis.text = element_text(size = 20), axis.title.x = element_text(size = 16, vjust = 2),
        axis.title.y = element_blank(), axis.text.y = element_blank(),
        legend.key.size = unit(1, "in"), legend.text = element_text(size = 18),
        legend.title = element_text(size = 20, hjust = 0.5), axis.ticks.length = unit(0.1, "in"),
        legend.position = "top", legend.title.position = "top") +
  coord_cartesian(ylim = c(0, 4.2), xlim = c(0,21))

# Create plots for Evergreen Needleleaf

abundPlot_Internal_EvergreenNeedleleaf <- ggplot() +
  scale_color_manual(values = plasma(3, begin = 0, end = 0.75, direction = -1)) +
  labs(x = "\n Evergreen Needleleaf Land Cover (%)", colour = "Region") +
  geom_point(data = internal, aes(x = pland_c01_evergreen_needleleaf, y = prediction, colour = loc), alpha = 0.01, size = 2, show.legend = FALSE) +
  geom_smooth(data = internal, aes(x = pland_c01_evergreen_needleleaf, y = prediction, group = boot_loc), colour = "grey", alpha = 0.5, method = "gam", se = T) +
  geom_smooth(data = internal, aes(x = pland_c01_evergreen_needleleaf, y = prediction, colour = loc), method = "gam", se = T, size = 3) +
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +
  theme(axis.text = element_text(size = 20), axis.title.x = element_text(size = 24, vjust = 2),
        axis.text.y = element_blank(), axis.title.y = element_blank(),
        legend.key.size = unit(1, "in"), legend.text = element_text(size = 20),
        legend.title = element_text(size = 20), axis.ticks.length = unit(0.1, "in")) +
  coord_cartesian(ylim = c(0, 4.2), xlim = c(0,105))

abundPlot_ExternalOut_EvergreenNeedleleaf <- ggplot() +
  scale_color_manual(values = plasma(3, begin = 0, end = 0.75, direction = -1)) +
  labs(x = "\n Evergreen Needleleaf Land Cover (%)", y = "Predicted CBCH Rel. Abund. \n", colour = "Scenario") +
  geom_point(data = external_out, aes(x = pland_c01_evergreen_needleleaf, y = prediction, colour = scenario), alpha = 0.01, size = 2, show.legend = FALSE) +
  geom_smooth(data = external_out, aes(x = pland_c01_evergreen_needleleaf, y = prediction, group = boot_loc), colour = "grey", alpha = 0.5, method = "gam", se = T) +
  geom_smooth(data = external_out, aes(x = pland_c01_evergreen_needleleaf, y = prediction, colour = scenario), method = "gam", se = T, size = 3) +
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +
  theme(axis.text = element_text(size = 20), axis.title.x = element_text(size = 16, vjust = 2),
        axis.title.y = element_text(size = 16),
        legend.key.size = unit(1, "in"), legend.text = element_text(size = 20),
        legend.title = element_text(size = 20), axis.ticks.length = unit(0.1, "in"),
        legend.position = "none") +
  coord_cartesian(ylim = c(0, 4.2), xlim = c(0,105))

abundPlot_ExternalIn_EvergreenNeedleleaf <- ggplot() +
  scale_color_manual(values = plasma(3, begin = 0, end = 0.75, direction = -1)) +
  labs(x = "\n Evergreen Needleleaf Land Cover (%)", colour = "Scenario") +
  geom_point(data = external_in, aes(x = pland_c01_evergreen_needleleaf, y = prediction, colour = scenario), alpha = 0.01, size = 2, show.legend = FALSE) +
  geom_smooth(data = external_in, aes(x = pland_c01_evergreen_needleleaf, y = prediction, group = boot_loc), colour = "grey", alpha = 0.5, method = "gam", se = T) +
  geom_smooth(data = external_in, aes(x = pland_c01_evergreen_needleleaf, y = prediction, colour = scenario), method = "gam", se = T, size = 3) +
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +
  theme(axis.text = element_text(size = 20), axis.title.x = element_text(size = 16, vjust = 2),
        axis.title.y = element_blank(), axis.text.y = element_blank(),
        legend.key.size = unit(1, "in"), legend.text = element_text(size = 20),
        legend.title = element_text(size = 20), axis.ticks.length = unit(0.1, "in"),
        legend.position = "none") +
  coord_cartesian(ylim = c(0, 4.2), xlim = c(0,105))

# Create title Grob to append to figures

plot_y_title_pred <- textGrob("Predicted CBCH Rel. Abund. \n", hjust = 0.4,
                         vjust = 0.8, rot = 90, gp = gpar(fontsize = 24))

# Gather internal plots, add y axis title, and save.

# Urban-related variables

plot_internal_urban <- ggarrange(abundPlot_Internal_RoadDensity,
                                 abundPlot_Internal_NighttimeLight,
                                 abundPlot_Internal_UrbanCover,
                                 nrow = 1,
                                 ncol = 3,
                                 common.legend = TRUE,
                                 legend = "top")

ggsave(plot_internal_urban, filename = "./Outputs/Plots/Figure2_Top.pdf",
         width = 7500, height = 3000, units = "px")

rm(plot_internal_urban)

# Forest-related variables

plot_internal_forest <- ggarrange(abundPlot_Internal_CanopyHeight,
                                  abundPlot_Internal_EvergreenNeedleleaf,
                                  nrow = 1,
                                  ncol = 2,
                                  common.legend = TRUE,
                                  legend = "top")

ggsave(plot_internal_forest, filename = "./Outputs/Plots/Figure2_Bottom.pdf",
       width = 5000, height = 3000, units = "px")

rm(plot_internal_forest)

# Gather external plots, unite figure legends, add y-axis label, and save.

# Urban-related Variables

# Grab legends

legend_out <- get_legend(abundPlot_ExternalOut_RoadDensity)
legend_in <- get_legend(abundPlot_ExternalIn_RoadDensity)

rm_legend <- function(p){p + theme(legend.position = "none")}

abundPlot_ExternalOut_RoadDensity <- rm_legend(abundPlot_ExternalOut_RoadDensity)
abundPlot_ExternalIn_RoadDensity <- rm_legend(abundPlot_ExternalIn_RoadDensity)

# Arrange plots and legends

plots_urb_top <- ggarrange(abundPlot_ExternalOut_RoadDensity,
                           abundPlot_ExternalIn_RoadDensity,
                           nrow = 1, ncol = 2)

legends <- ggarrange(legend_out, legend_in, nrow = 1)

# Unite plots and legends and save.

plot_external_urban_top <- ggarrange(legends, plots_urb_top, heights = c(0.25, 0.75), nrow = 2)

ggsave(plot_external_urban_top, filename = "./Outputs/Plots/Figure5_Top.pdf",
       width = 5000, height = 4000, units = "px")

rm(plots_urb_top)
rm(plot_external_urban_top)

plots_urb_middle <- ggarrange(abundPlot_ExternalOut_NighttimeLight,
                              abundPlot_ExternalIn_NighttimeLight,
                              nrow = 1, ncol = 2)

ggsave(plots_urb_middle, filename = "./Outputs/Plots/Figure5_Middle.pdf",
         width = 5000, height = 3000, units = "px")

rm(plots_urb_middle)

plots_urb_bottom <- ggarrange(abundPlot_ExternalOut_UrbanCover,
                              abundPlot_ExternalIn_UrbanCover,
                              nrow = 1, ncol = 2)

ggsave(plots_urb_bottom, filename = "./Outputs/Plots/Figure5_Bottom.pdf",
       width = 5000, height = 3000, units = "px")

rm(plots_urb_bottom)

# Forest-related variables

# Arrange plots

abundPlot_ExternalOut_CanopyHeight <- rm_legend(abundPlot_ExternalOut_CanopyHeight)
abundPlot_ExternalIn_CanopyHeight <- rm_legend(abundPlot_ExternalIn_CanopyHeight)


plots_forest_top <- ggarrange(abundPlot_ExternalOut_CanopyHeight,
                             abundPlot_ExternalIn_CanopyHeight,
                             ncol = 2, nrow = 1)

# Unite plots and legends and save.

plot_external_forest_top <- ggarrange(legends, plots_forest_top, heights = c(0.25, 0.75), nrow = 2)

ggsave(plot_external_forest_top, filename = "./Outputs/Plots/Figure6_Top.pdf",
       width = 5000, height = 4000, units = "px")

rm(plots_forest_top)
rm(plot_external_forest_top)

plots_forest_bottom <- ggarrange(abundPlot_ExternalOut_EvergreenNeedleleaf,
                                 abundPlot_ExternalIn_EvergreenNeedleleaf,
                                 nrow = 1, ncol = 2)

ggsave(plots_forest_bottom, filename = "./Outputs/Plots/Figure6_Bottom.pdf",
       width = 5000, height = 3000, units = "px")

rm(plots_forest_bottom)

# Counties

abund_long_counties <- pivot_longer(abundance_predictions_counties, cols = starts_with("boot"), names_to = "bootstrap", values_to = "prediction")
abund_long_counties$boot_loc <- paste0(abund_long_counties$bootstrap, abund_long_counties$loc)
abund_long_counties$boot_loc_InOut <- paste0(abund_long_counties$bootstrap, abund_long_counties$loc_InOut)

# Add a column identifying which urban centre a cell is in.

abund_long_counties <- mutate(abund_long_counties, region = case_when(cell_id %in% wa_grid$cell_id ~ "WA",
                                                 cell_id %in% or_grid$cell_id ~ "OR",
                                                 cell_id %in% mv_grid$cell_id ~ "MV",
                                                 cell_id %in% vi_grid$cell_id ~ "VI",
                                                 cell_id %in% sf_grid$cell_id ~ "SF",
                                                 cell_id %in% nc_grid$cell_id ~ "NC"))

# Subset to create sets we'll use for plots

# Create 'internal' for use in Figure 3 - demonstrating the relationship to focal
# variables within the area a model was built on.

internal_counties <- filter(abund_long_counties, loc == region)

# Create 'external_out' for use in Figure 5 - comparing the real-world predicted
# abundance and the predicted abundance under imagined scenarios in allopatric areas 

external_counties_out <- abund_long_counties %>%
  filter(loc %in% as.factor("WA") & region %in% c("VI", "SF", "NC") | loc %in% as.factor("OR") %in% region %in% c("VI", "SF", "NC") |
           loc %in% as.factor("MV") & region %in% c("VI", "SF", "NC") | loc %in% as.factor("VI") & region == "VI" | loc %in% as.factor("SF") & region == "SF" |
           loc %in% as.factor("NC") & region == "NC") %>%
  mutate(scenario = case_when(loc %in% c("VI", "SF", "NC") ~ "Real-world",
                              loc %in% c("WA", "OR", "MV") ~ "BCCH Present"))

external_counties_out$scenario = factor(external_counties_out$scenario, c("Real-world", "BCCH Present"))

# Create 'external_in' for use in Figure 5 - comparing the real-world predicted
# abundance and the predicted abundance under imagined scenarios in sympatric areas

external_counties_in <- abund_long_counties %>%
  filter(loc %in% as.factor("VI") & region %in% c("WA", "OR", "MV") | loc %in% as.factor("SF") %in% region %in% c("WA", "OR", "MV") |
           loc %in% as.factor("NC") & region %in% c("WA", "OR", "MV") | loc %in% as.factor("WA") & region == "WA" | loc %in% as.factor("OR") & region == "OR" |
           loc %in% as.factor("MV") & region == "MV") %>%
  mutate(scenario = case_when(loc %in% c("WA", "OR", "MV") ~ "Real-world",
                              loc %in% c("VI", "SF", "NC") ~ "BCCH Absent"))

external_counties_in$scenario = factor(external_counties_in$scenario, c("Real-world", "BCCH Absent"))

counties_palette <- c(plasma(3, begin = 0.5, end = 1, direction = -1),
                     plasma(3, begin = 0, end = 0.5, direction = -1))

# Internal plots

# Road Density

abundPlot_Counties_Internal_RoadDensity <- ggplot() +
  scale_color_manual(values = plasma(3, begin = 0, end = 0.75, direction = -1)) +
  labs(x = "\n Road Density (m/km^2)",  y = "Predicted CBCH Rel. Abund. \n", colour = "Location") +
  geom_point(data = internal_counties, aes(x = road_dens, y = prediction, colour = loc_InOut), alpha = 0.01, size = 2, show.legend = FALSE) +
  geom_smooth(data = internal_counties, aes(x = road_dens, y = prediction, group = boot_loc_InOut), colour = "grey", alpha = 0.5, method = "gam", se = T) +
  geom_smooth(data = internal_counties, aes(x = road_dens, y = prediction, colour = loc_InOut), method = "gam", se = T, size = 3) +
  scale_x_continuous(expand = c(0, 0), breaks = c(5000, 10000, 15000)) + 
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +
  theme(axis.text = element_text(size = 20), axis.title.x = element_text(size = 24, vjust = 2),
        axis.title.y = element_text(size = 24),
        legend.key.size = unit(1, "in"), legend.text = element_text(size = 20),
        legend.title = element_text(size = 20), axis.ticks.length = unit(0.1, "in")) +
  coord_cartesian(ylim = c(0, 4.2), xlim = c(0,16000))

abundPlot_Counties_ExternalOut_RoadDensity <- ggplot() +
  scale_color_manual(values = plasma(3, begin = 0, end = 0.75, direction = -1)) +
  labs(x = "\n Road Density (m/km^2)", y = "Predicted CBCH Rel. Abund. \n", colour = "Allopatric Areas") +
  geom_point(data = external_counties_out, aes(x = road_dens, y = prediction, colour = scenario), alpha = 0.01, size = 2, show.legend = FALSE) +
  geom_smooth(data = external_counties_out, aes(x = road_dens, y = prediction, group = boot_loc_InOut), colour = "grey", alpha = 0.5, method = "gam", se = T) +
  geom_smooth(data = external_counties_out, aes(x = road_dens, y = prediction, colour = scenario), method = "gam", se = T, size = 3) +
  scale_x_continuous(expand = c(0, 0), breaks = c(5000, 10000, 15000)) + 
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +
  theme(axis.text = element_text(size = 20), axis.title.x = element_text(size = 16, vjust = 2),
        axis.title.y = element_text(size = 16),
        legend.key.size = unit(1, "in"), legend.text = element_text(size = 18),
        legend.title = element_text(size = 20, hjust = 0.5), axis.ticks.length = unit(0.1, "in"),
        legend.position = "top", legend.title.position = "top") +
  coord_cartesian(ylim = c(0, 4.2), xlim = c(0,16000))

abundPlot_Counties_ExternalIn_RoadDensity <- ggplot() +
  scale_color_manual(values = plasma(3, begin = 0, end = 0.75, direction = -1)) +
  labs(x = "\n Road Density (m/km^2)", colour = "Sympatric Areas") +
  geom_point(data = external_counties_in, aes(x = road_dens, y = prediction, colour = scenario), alpha = 0.01, size = 2, show.legend = FALSE) +
  geom_smooth(data = external_counties_in, aes(x = road_dens, y = prediction, group = boot_loc_InOut), colour = "grey", alpha = 0.5, method = "gam", se = T) +
  geom_smooth(data = external_counties_in, aes(x = road_dens, y = prediction, colour = scenario), method = "gam", se = T, size = 3) +
  scale_x_continuous(expand = c(0, 0), breaks = c(5000, 10000, 15000)) + 
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +
  theme(axis.text = element_text(size = 20), axis.title.x = element_text(size = 16, vjust = 2),
        axis.title.y = element_blank(), axis.text.y = element_blank(),
        legend.key.size = unit(1, "in"), legend.text = element_text(size = 18),
        legend.title = element_text(size = 20, hjust = 0.5), axis.ticks.length = unit(0.1, "in"),
        legend.position = "top", legend.title.position = "top") +
  coord_cartesian(ylim = c(0, 4.2), xlim = c(0,16000))

# Nighttime Light

abundPlot_Counties_Internal_NighttimeLight <- ggplot() +
  scale_color_manual(values = plasma(3, begin = 0, end = 0.75, direction = -1)) +
  labs(x = "\n Nighttime Light Radiance", colour = "Location") +
  geom_point(data = internal_counties, aes(x = nighttime_light, y = prediction, colour = loc_InOut), alpha = 0.01, size = 2, show.legend = FALSE) +
  geom_smooth(data = internal_counties, aes(x = nighttime_light, y = prediction, group = boot_loc_InOut), colour = "grey", alpha = 0.5, method = "gam", se = T) +
  geom_smooth(data = internal_counties, aes(x = nighttime_light, y = prediction, colour = loc_InOut), method = "gam", se = T, size = 3) +
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +
  theme(axis.text = element_text(size = 20), axis.title.x = element_text(size = 24, vjust = 2),
        axis.text.y = element_blank(), axis.title.y = element_blank(),
        legend.key.size = unit(1, "in"), legend.text = element_text(size = 20),
        legend.title = element_text(size = 20), axis.ticks.length = unit(0.1, "in")) +
  coord_cartesian(ylim = c(0, 4.2), xlim = c(0,35))

abundPlot_Counties_ExternalOut_NighttimeLight <- ggplot() +
  scale_color_manual(values = plasma(3, begin = 0, end = 0.75, direction = -1)) +
  labs(x = "\n Nighttime Light Radiance",  y = "Predicted CBCH Rel. Abund. \n", colour = "Location") +
  geom_point(data = external_counties_out, aes(x = nighttime_light, y = prediction, colour = scenario), alpha = 0.01, size = 2, show.legend = FALSE) +
  geom_smooth(data = external_counties_out, aes(x = nighttime_light, y = prediction, group = boot_loc_InOut), colour = "grey", alpha = 0.5, method = "gam", se = T) +
  geom_smooth(data = external_counties_out, aes(x = nighttime_light, y = prediction, colour = scenario), method = "gam", se = T, size = 3) +
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +
  theme(axis.text = element_text(size = 20), axis.title.x = element_text(size = 16, vjust = 2),
        axis.title.y = element_text(size = 16),
        legend.key.size = unit(1, "in"), legend.text = element_text(size = 20),
        legend.title = element_text(size = 20), axis.ticks.length = unit(0.1, "in"),
        legend.position = "none") +
  coord_cartesian(ylim = c(0, 4.2), xlim = c(0,35))

abundPlot_Counties_ExternalIn_NighttimeLight <- ggplot() +
  scale_color_manual(values = plasma(3, begin = 0, end = 0.75, direction = -1)) +
  labs(x = "\n Nighttime Light Radiance", colour = "Location") +
  geom_point(data = external_counties_in, aes(x = nighttime_light, y = prediction, colour = scenario), alpha = 0.01, size = 2, show.legend = FALSE) +
  geom_smooth(data = external_counties_in, aes(x = nighttime_light, y = prediction, group = boot_loc_InOut), colour = "grey", alpha = 0.5, method = "gam", se = T) +
  geom_smooth(data = external_counties_in, aes(x = nighttime_light, y = prediction, colour = scenario), method = "gam", se = T, size = 3) +
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +
  theme(axis.text = element_text(size = 20), axis.title.x = element_text(size = 16, vjust = 2),
        axis.title.y = element_blank(), axis.text.y = element_blank(),
        legend.key.size = unit(1, "in"), legend.text = element_text(size = 20),
        legend.title = element_text(size = 20), axis.ticks.length = unit(0.1, "in"),
        legend.position = "none") +
  coord_cartesian(ylim = c(0, 4.2), xlim = c(0,35))

# Urban Cover

abundPlot_Counties_Internal_UrbanCover <- ggplot() +
  scale_color_manual(values = plasma(3, begin = 0, end = 0.75, direction = -1)) +
  labs(x = "\n Urban Land Cover (%)", colour = "Location") +
  geom_point(data = internal_counties, aes(x = pland_c13_urban, y = prediction, colour = loc_InOut), alpha = 0.01, size = 2, show.legend = FALSE) +
  geom_smooth(data = internal_counties, aes(x = pland_c13_urban, y = prediction, group = boot_loc_InOut), colour = "grey", alpha = 0.5, method = "gam", se = T) +
  geom_smooth(data = internal_counties, aes(x = pland_c13_urban, y = prediction, colour = loc_InOut), method = "gam", se = T, size = 3) +
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +
  theme(axis.text = element_text(size = 20), axis.title.x = element_text(size = 24, vjust = 2),
        axis.text.y = element_blank(), axis.title.y = element_blank(),
        legend.key.size = unit(1, "in"), legend.text = element_text(size = 20),
        legend.title = element_text(size = 20), axis.ticks.length = unit(0.1, "in")) +
  coord_cartesian(ylim = c(0, 4.2), xlim = c(0,105))

abundPlot_Counties_ExternalOut_UrbanCover <- ggplot() +
  scale_color_manual(values = plasma(3, begin = 0, end = 0.75, direction = -1)) +
  labs(x = "\n Urban Land Cover (%)",  y = "Predicted CBCH Rel. Abund. \n", colour = "Location") +
  geom_point(data = external_counties_out, aes(x = pland_c13_urban, y = prediction, colour = scenario), alpha = 0.01, size = 2, show.legend = FALSE) +
  geom_smooth(data = external_counties_out, aes(x = pland_c13_urban, y = prediction, group = boot_loc_InOut), colour = "grey", alpha = 0.5, method = "gam", se = T) +
  geom_smooth(data = external_counties_out, aes(x = pland_c13_urban, y = prediction, colour = scenario), method = "gam", se = T, size = 3) +
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +
  theme(axis.text = element_text(size = 20), axis.title.x = element_text(size = 16, vjust = 2),
        axis.title.y = element_text(size = 16),
        legend.key.size = unit(1, "in"), legend.text = element_text(size = 20),
        legend.title = element_text(size = 20), axis.ticks.length = unit(0.1, "in"),
        legend.position = "none") +
  coord_cartesian(ylim = c(0, 4.2), xlim = c(0,105))

abundPlot_Counties_ExternalIn_UrbanCover <- ggplot() +
  scale_color_manual(values = plasma(3, begin = 0, end = 0.75, direction = -1)) +
  labs(x = "\n Urban Land Cover (%)", colour = "Location") +
  geom_point(data = external_counties_in, aes(x = pland_c13_urban, y = prediction, colour = scenario), alpha = 0.01, size = 2, show.legend = FALSE) +
  geom_smooth(data = external_counties_in, aes(x = pland_c13_urban, y = prediction, group = boot_loc_InOut), colour = "grey", alpha = 0.5, method = "gam", se = T) +
  geom_smooth(data = external_counties_in, aes(x = pland_c13_urban, y = prediction, colour = scenario), method = "gam", se = T, size = 3) +
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +
  theme(axis.text = element_text(size = 20), axis.title.x = element_text(size = 16, vjust = 2),
        axis.title.y = element_blank(), axis.text.y = element_blank(),
        legend.key.size = unit(1, "in"), legend.text = element_text(size = 20),
        legend.title = element_text(size = 20), axis.ticks.length = unit(0.1, "in"),
        legend.position = "none") +
  coord_cartesian(ylim = c(0, 4.2), xlim = c(0,105))


# Canopy Height

abundPlot_Counties_Internal_CanopyHeight <- ggplot() +
  scale_color_manual(values = plasma(3, begin = 0, end = 0.75, direction = -1)) +
  labs(x = "\n Mean Canopy Height (m)",  y = "Predicted CBCH Rel. Abund. \n", colour = "Location") +
  geom_point(data = internal_counties, aes(x = canopy_height, y = prediction, colour = loc_InOut), alpha = 0.01, size = 2, show.legend = FALSE) +
  geom_smooth(data = internal_counties, aes(x = canopy_height, y = prediction, group = boot_loc_InOut), colour = "grey", alpha = 0.5, method = "gam", se = T) +
  geom_smooth(data = internal_counties, aes(x = canopy_height, y = prediction, colour = loc_InOut), method = "gam", se = T, size = 3) +
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +
  theme(axis.text = element_text(size = 20), axis.title.x = element_text(size = 24, vjust = 2),
        axis.title.y = element_text(size = 24),
        legend.key.size = unit(1, "in"), legend.text = element_text(size = 20),
        legend.title = element_text(size = 20), axis.ticks.length = unit(0.1, "in")) +
  coord_cartesian(ylim = c(0, 4.2), xlim = c(0,21))

abundPlot_Counties_ExternalOut_CanopyHeight <- ggplot() +
  scale_color_manual(values = plasma(3, begin = 0, end = 0.75, direction = -1)) +
  labs(x = "\n Mean Canopy Height (m)",  y = "Predicted CBCH Rel. Abund. \n", colour = "Allopatric Areas") +
  geom_point(data = external_counties_out, aes(x = canopy_height, y = prediction, colour = scenario), alpha = 0.01, size = 2, show.legend = FALSE) +
  geom_smooth(data = external_counties_out, aes(x = canopy_height, y = prediction, group = boot_loc_InOut), colour = "grey", alpha = 0.5, method = "gam", se = T) +
  geom_smooth(data = external_counties_out, aes(x = canopy_height, y = prediction, colour = scenario), method = "gam", se = T, size = 3) +
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +
  theme(axis.text = element_text(size = 20), axis.title.x = element_text(size = 16, vjust = 2),
        axis.title.y = element_text(size = 16),
        legend.key.size = unit(1, "in"), legend.text = element_text(size = 18),
        legend.title = element_text(size = 20, hjust = 0.5), axis.ticks.length = unit(0.1, "in"),
        legend.position = "top", legend.title.position = "top") +
  coord_cartesian(ylim = c(0, 4.2), xlim = c(0,21))

abundPlot_Counties_ExternalIn_CanopyHeight <- ggplot() +
  scale_color_manual(values = plasma(3, begin = 0, end = 0.75, direction = -1)) +
  labs(x = "\n Mean Canopy Height (m)", colour = "Sympatric Areas") +
  geom_point(data = external_counties_in, aes(x = canopy_height, y = prediction, colour = scenario), alpha = 0.01, size = 2, show.legend = FALSE) +
  geom_smooth(data = external_counties_in, aes(x = canopy_height, y = prediction, group = boot_loc_InOut), colour = "grey", alpha = 0.5, method = "gam", se = T) +
  geom_smooth(data = external_counties_in, aes(x = canopy_height, y = prediction, colour = scenario), method = "gam", se = T, size = 3) +
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +
  theme(axis.text = element_text(size = 20), axis.title.x = element_text(size = 16, vjust = 2),
        axis.title.y = element_blank(), axis.text.y = element_blank(),
        legend.key.size = unit(1, "in"), legend.text = element_text(size = 18),
        legend.title = element_text(size = 20, hjust = 0.5), axis.ticks.length = unit(0.1, "in"),
        legend.position = "top", legend.title.position = "top") +
  coord_cartesian(ylim = c(0, 4.2), xlim = c(0,21))

# Evergreen Needleleaf

abundPlot_Counties_Internal_EvergreenNeedleleaf <- ggplot() +
  scale_color_manual(values = plasma(3, begin = 0, end = 0.75, direction = -1)) +
  labs(x = "\n Evergreen Needleleaf Land Cover (%)", colour = "Location") +
  geom_point(data = internal_counties, aes(x = pland_c01_evergreen_needleleaf, y = prediction, colour = loc_InOut), alpha = 0.01, size = 2, show.legend = FALSE) +
  geom_smooth(data = internal_counties, aes(x = pland_c01_evergreen_needleleaf, y = prediction, group = boot_loc_InOut), colour = "grey", alpha = 0.5, method = "gam", se = T) +
  geom_smooth(data = internal_counties, aes(x = pland_c01_evergreen_needleleaf, y = prediction, colour = loc_InOut), method = "gam", se = T, size = 3) +
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +
  theme(axis.text = element_text(size = 20), axis.title.x = element_text(size = 24, vjust = 2),
        axis.text.y = element_blank(), axis.title.y = element_blank(),
        legend.key.size = unit(1, "in"), legend.text = element_text(size = 20),
        legend.title = element_text(size = 20), axis.ticks.length = unit(0.1, "in")) +
  coord_cartesian(ylim = c(0, 4.2), xlim = c(0,105))

abundPlot_Counties_ExternalOut_EvergreenNeedleleaf <- ggplot() +
  scale_color_manual(values = plasma(3, begin = 0, end = 0.75, direction = -1)) +
  labs(x = "\n Evergreen Needleleaf Land Cover (%)", y = "Predicted CBCH Rel. Abund. \n", colour = "Scenario") +
  geom_point(data = external_counties_out, aes(x = pland_c01_evergreen_needleleaf, y = prediction, colour = scenario), alpha = 0.01, size = 2, show.legend = FALSE) +
  geom_smooth(data = external_counties_out, aes(x = pland_c01_evergreen_needleleaf, y = prediction, group = boot_loc_InOut), colour = "grey", alpha = 0.5, method = "gam", se = T) +
  geom_smooth(data = external_counties_out, aes(x = pland_c01_evergreen_needleleaf, y = prediction, colour = scenario), method = "gam", se = T, size = 3) +
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +
  theme(axis.text = element_text(size = 20), axis.title.x = element_text(size = 16, vjust = 2),
        axis.title.y = element_text(size = 16),
        legend.key.size = unit(1, "in"), legend.text = element_text(size = 20),
        legend.title = element_text(size = 20), axis.ticks.length = unit(0.1, "in"),
        legend.position = "none") +
  coord_cartesian(ylim = c(0, 4.2), xlim = c(0,105))

abundPlot_Counties_ExternalIn_EvergreenNeedleleaf <- ggplot() +
  scale_color_manual(values = plasma(3, begin = 0, end = 0.75, direction = -1)) +
  labs(x = "\n Evergreen Needleleaf Land Cover (%)", colour = "Scenario") +
  geom_point(data = external_counties_in, aes(x = pland_c01_evergreen_needleleaf, y = prediction, colour = scenario), alpha = 0.01, size = 2, show.legend = FALSE) +
  geom_smooth(data = external_counties_in, aes(x = pland_c01_evergreen_needleleaf, y = prediction, group = boot_loc_InOut), colour = "grey", alpha = 0.5, method = "gam", se = T) +
  geom_smooth(data = external_counties_in, aes(x = pland_c01_evergreen_needleleaf, y = prediction, colour = scenario), method = "gam", se = T, size = 3) +
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +
  theme(axis.text = element_text(size = 20), axis.title.x = element_text(size = 16, vjust = 2),
        axis.title.y = element_blank(), axis.text.y = element_blank(),
        legend.key.size = unit(1, "in"), legend.text = element_text(size = 20),
        legend.title = element_text(size = 20), axis.ticks.length = unit(0.1, "in"),
        legend.position = "none") +
  coord_cartesian(ylim = c(0, 4.2), xlim = c(0,105))

# Gather internal plots, add y axis title, and save.

# Urban-related variables

plot_internal_urban_counties <- ggarrange(abundPlot_Counties_Internal_RoadDensity,
                                          abundPlot_Counties_Internal_NighttimeLight,
                                          abundPlot_Counties_Internal_UrbanCover,
                                          nrow = 1,
                                          ncol = 3,
                                          common.legend = TRUE,
                                          legend = "top")

ggsave(plot_internal_urban_counties, filename = "./Outputs/Plots/Figure1_Top.pdf",
         width = 7500, height = 3000, units = "px")

rm(plot_internal_urban_counties)

# Forest-related variables

plot_internal_forest_counties <- ggarrange(abundPlot_Counties_Internal_CanopyHeight,
                                  abundPlot_Counties_Internal_EvergreenNeedleleaf,
                                  nrow = 1,
                                  ncol = 2,
                                  common.legend = TRUE,
                                  legend = "top")

ggsave(plot_internal_forest_counties, filename = "./Outputs/Plots/Figure1_Bottom.pdf",
         width = 5000, height = 3000, units = "px")

rm(plot_internal_forest_counties)

# Gather external plots, unite figure legends, add y-axis label, and save.

# Urban-related Variables

# Grab legends

legend_out <- get_legend(abundPlot_Counties_ExternalOut_RoadDensity)
legend_in <- get_legend(abundPlot_Counties_ExternalIn_RoadDensity)

rm_legend <- function(p){p + theme(legend.position = "none")}

abundPlot_Counties_ExternalOut_RoadDensity <- rm_legend(abundPlot_Counties_ExternalOut_RoadDensity)
abundPlot_Counties_ExternalIn_RoadDensity <- rm_legend(abundPlot_Counties_ExternalIn_RoadDensity)

# Arrange plots and legends

# Top panel

plots_urb_counties_top <- ggarrange(abundPlot_Counties_ExternalOut_RoadDensity,
                                abundPlot_Counties_ExternalIn_RoadDensity,
                                nrow = 1, ncol = 2)

legends <- ggarrange(legend_out, legend_in, nrow = 1)

plot_external_urban_counties_top <- ggarrange(legends, plots_urb_counties_top, heights = c(0.25, 0.75), nrow = 2)


ggsave(plot_external_urban_counties_top, filename = "./Outputs/Plots/Figure3_Top.pdf",
       width = 5000, height = 4000, units = "px")

rm(plots_urb_counties_top)
rm(plot_external_urban_counties_top)


# Middle panel

plots_urb_counties_middle <- ggarrange(abundPlot_Counties_ExternalOut_NighttimeLight,
                                       abundPlot_Counties_ExternalIn_NighttimeLight,
                                       nrow = 1, ncol = 2)

ggsave(plots_urb_counties_middle, filename = "./Outputs/Plots/Figure3_Middle.pdf",
       width = 5000, height = 3000, units = "px")

rm(plots_urb_counties_middle)

# Bottom panel

plots_urb_counties_bottom <- ggarrange(abundPlot_Counties_ExternalOut_UrbanCover,
                                       abundPlot_Counties_ExternalIn_UrbanCover,
                                       nrow = 1, ncol = 2)

ggsave(plots_urb_counties_bottom, filename = "./Outputs/Plots/Figure3_Bottom.pdf",
       width = 5000, height = 3000, units = "px")

rm(plots_urb_counties_bottom)


# Forest-related variables

# Arrange plots

abundPlot_Counties_ExternalOut_CanopyHeight <- rm_legend(abundPlot_Counties_ExternalOut_CanopyHeight)
abundPlot_Counties_ExternalIn_CanopyHeight <- rm_legend(abundPlot_Counties_ExternalIn_CanopyHeight)

plots_forest_counties_top <- ggarrange(abundPlot_Counties_ExternalOut_CanopyHeight,
                                       abundPlot_Counties_ExternalIn_CanopyHeight,
                                       nrow = 1, ncol = 2)


# Unite plots and legends and save.

plot_external_forest_counties_top <- ggarrange(legends, plots_forest_counties_top, heights = c(0.25, 0.75), nrow = 2)

ggsave(plot_external_forest_counties_top, filename = "./Outputs/Plots/Figure4_Top.pdf",
         width = 5000, height = 4000, units = "px")

rm(plots_forest_counties_top)
rm(plot_external_forest_counties_top)
                         
plots_forest_counties_bottom <- ggarrange(abundPlot_Counties_ExternalOut_EvergreenNeedleleaf,
                                          abundPlot_Counties_ExternalIn_EvergreenNeedleleaf,
                                          nrow = 1, ncol = 2)

ggsave(plots_forest_counties_bottom, filename = "./Outputs/Plots/Figure4_Bottom.pdf",
       width = 5000, height = 3000, units = "px")

rm(plots_forest_counties_bottom)

####################### PREDICTIVE PERFORMANCE METRICS #########################

ppms <- data.frame(metric = rep(c("Kappa", "PR_AUC", "ROC_AUC", "F1", "Count_Spearman"), each = length(locs)),
                   loc = rep(c("Allopatric", "Sympatric", "Region-wide", "SEA", "POR", "MV", "SVI", "SF", "NC"), times = 5))

for(i in seq_along(seeds)) {
  
  for(j in ppms$metric) {
    
    out_metric <- out_assessment[[i]][j]
    ppms[ppms$metric == j & ppms$loc == "Allopatric", paste0("rank_boot", seeds[i])] <- out_metric
    
    in_metric <- in_assessment[[i]][j]
    ppms[ppms$metric == j & ppms$loc == "Sympatric", paste0("rank_boot", seeds[i])] <- in_metric
    
    all_metric <- all_assessment[[i]][j]
    ppms[ppms$metric == j & ppms$loc == "Region-wide", paste0("rank_boot", seeds[i])] <- all_metric
    
    wa_metric <- counties_assessment$wa[[i]][[j]]
    ppms[ppms$metric == j & ppms$loc == "SEA", paste0("rank_boot", seeds[i])] <- wa_metric
    
    or_metric <- counties_assessment$or[[i]][[j]]
    ppms[ppms$metric == j & ppms$loc == "POR", paste0("rank_boot", seeds[i])] <- or_metric
    
    mv_metric <- counties_assessment$mv[[i]][[j]]
    ppms[ppms$metric == j & ppms$loc == "MV", paste0("rank_boot", seeds[i])] <- mv_metric
    
    vi_metric <- counties_assessment$vi[[i]][[j]]
    ppms[ppms$metric == j & ppms$loc == "SVI", paste0("rank_boot", seeds[i])] <- vi_metric
    
    sf_metric <- counties_assessment$sf[[i]][[j]]
    ppms[ppms$metric == j & ppms$loc == "SF", paste0("rank_boot", seeds[i])] <- sf_metric
    
    nc_metric <- counties_assessment$nc[[i]][[j]]
    ppms[ppms$metric == j & ppms$loc == "NC", paste0("rank_boot", seeds[i])] <- nc_metric
    
  }
  
}

ppms$metric[ppms$metric == "Count_Spearman"] <- "SRC"

ppms$metric <- factor(ppms$metric, levels = c("ROC_AUC", "PR_AUC", "Kappa",
                                              "F1", "SRC"))

ppms$loc <- factor(ppms$loc, levels = c("Allopatric", "Sympatric", "Region-wide",
                                        "SVI", "SF", "NC", "MV", "SEA", "POR"))

ppm_plot <- ppms[ppms$loc %in% c("Sympatric", "Allopatric"),] %>%
  pivot_longer(3:length(.), names_to = "boot", values_to = "value") %>%
  ggplot(aes(x = metric, y = value, fill = loc)) +
  geom_violin(position = "dodge") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  labs(fill = "Location") + 
  xlab("") +
  ylab("") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  facet_grid(. ~ metric, scales = "free", space = "free")

ppm_counties_plot <- ppms[ppms$loc %in% c("SEA", "POR", "MV", "SVI", "SF", "NC"),] %>%
  pivot_longer(3:length(.), names_to = "boot", values_to = "value") %>%
  ggplot(aes(x = metric, y = value, fill = loc)) +
  geom_violin(position = "dodge") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  labs(fill = "Location") + 
  xlab("") +
  ylab("") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  facet_grid(. ~ metric, scales = "free", space = "free")

ggsave("./Outputs/Plots/PPMs.png",
       ppm_plot,
       device = "png",
       width = 7,
       height = 5,
       units = "in")

ggsave("./Outputs/Plots/CountiesPPMs.png",
       ppm_counties_plot,
       device = "png",
       width = 7,
       height = 5,
       units = "in")
