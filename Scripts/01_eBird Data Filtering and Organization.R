################################################################################

# Script Title: Filtering and Visualizing the eBird Basic Dataset for Bird
# Conservation Region 5 wide analysis of the impact of Black-capped Chickadee
# Presence on Chestnut-backed Chickadee Relative Abundance at two spatial
# grains.

# Note: Much of the code in this script is taken from Dr. Matthew Strimas-Mackey's
# fantastic guide on the analysis of eBird Data, which can be found here:
# https://cornelllabofornithology.github.io/ebird-best-practices/.

################################################################################

## Install required packages (if not already installed).

# install.packages("librarian")

librarian::shelf(auk, lubridate, sf, gridExtra, tidyverse, utils)

## Resolve namespace conflicts

select <- dplyr::select

# In agreement with the eBird Dataset Terms of Use, original data files as
# downloaded from the eBird portal are not included in this repo. See
# https://ebird.org/data/download to make a request. We downloaded the
# eBird Basic Dataset Sampling Event Data and eBird Basic Dataset filtered
# to Chestnut-backed Chickadees (Poecile rufescens) in 5 files for Alaska,
# British Columbia, Washington, Oregon, and California respectively from the
# Apr-2023 release of the eBird Dataset & the same 5 files for Black-capped 
# Chickadee (Poecile atricapillus). We filtered data to Jan 2000 and later.
# See README for full citation.

### Filtering Chestnut-backed Chickadee Data

## Load downloaded eBird Basic Dataset for Washington and sampling database.

CBCH_WA <- auk_ebd("./Data/Raw/eBird-Basic-Dataset/CBCH/ebd_US-WA_chbchi_relApr-2023.txt", 
                   sep = "\t",
                   file_sampling = "./Data/Raw/EBD_Sampling/ebd_sampling_relApr-2023.txt")

# Filter to BCR 5, complete checklists, and Stationary and Traveling methods.

CBCH_filters_WA <- CBCH_WA %>%
  auk_bcr(bcr = c(5,32)) %>%
  auk_protocol(protocol = c("Stationary", "Traveling")) %>%
  auk_complete()

auk_filter(CBCH_filters_WA, 
           file = "./Data/Cleaned/eBird-Basic-Dataset/Filtered/CBCH/WA_CBCH_f_relApr-2023.txt", 
           file_sampling = "./Data/Raw/EBD_Sampling/sampling_f_relApr-2023.txt",
           overwrite = TRUE)

# Load filtered file back into R environment.

CBCH_WA_f <- read_ebd("./Data/Cleaned/eBird-Basic-Dataset/Filtered/CBCH/WA_CBCH_f_relApr-2023.txt",
                      sep = "\t")

## Repeat what was applied to Washington with Oregon, using the filtered sampling
## dataset we created above.

CBCH_OR <- auk_ebd("./Data/Raw/eBird-Basic-Dataset/CBCH/ebd_US-OR_chbchi_relApr-2023.txt", sep = "\t")

CBCH_filters_OR <- CBCH_OR %>%
  auk_bcr(bcr = 5) %>%
  auk_protocol(protocol = c("Stationary", "Traveling")) %>%
  auk_complete()

auk_filter(CBCH_filters_OR, file = "./Data/Cleaned/eBird-Basic-Dataset/Filtered/CBCH/OR_CBCH_f_relApr-2023.txt",
           overwrite = TRUE)

CBCH_OR_f <- read_ebd("./Data/Cleaned/eBird-Basic-Dataset/Filtered/CBCH/OR_CBCH_f_relApr-2023.txt", 
                      sep = "\t")

## Repeat with British Columbia.

CBCH_BC <- auk_ebd("./Data/Raw/eBird-Basic-Dataset/CBCH/ebd_CA-BC_chbchi_relApr-2023.txt", sep = "\t")

CBCH_filters_BC <- CBCH_BC %>%
  auk_bcr(bcr = 5) %>%
  auk_protocol(protocol = c("Stationary", "Traveling")) %>%
  auk_complete()

auk_filter(CBCH_filters_BC, file = "./Data/Cleaned/eBird-Basic-Dataset/Filtered/CBCH/BC_CBCH_f_relApr-2023.txt",
           overwrite = TRUE)

CBCH_BC_f <- read_ebd("./Data/Cleaned/eBird-Basic-Dataset/Filtered/CBCH/BC_CBCH_f_relApr-2023.txt", 
                      sep = "\t")

## Repeat with California.

CBCH_CA <- auk_ebd("./Data/Raw/eBird-Basic-Dataset/CBCH/ebd_US-CA_chbchi_relApr-2023.txt", sep = "\t")

CBCH_filters_CA <- CBCH_CA %>%
  auk_bcr(bcr = c(5, 32)) %>%
  auk_protocol(protocol = c("Stationary", "Traveling")) %>%
  auk_complete()

auk_filter(CBCH_filters_CA, file = "./Data/Cleaned/eBird-Basic-Dataset/Filtered/CBCH/CA_CBCH_f_relApr-2023.txt",
           overwrite = TRUE)

CBCH_CA_f <- read_ebd("./Data/Cleaned/eBird-Basic-Dataset/Filtered/CBCH/CA_CBCH_f_relApr-2023.txt", 
                      sep = "\t")

CBCH_CA_f <- CBCH_CA_f %>%
  filter(county %in% c(unique(CBCH_CA_f$county[CBCH_CA_f$bcr_code == 5]), 
                       "Napa", "Marin", "Solano", "Alameda", "San Mateo", "San Francisco",
                       "Contra Costa", "Santa Clara"))

## Compile into one object

CBCH_f <- rbind(CBCH_BC_f, CBCH_WA_f, CBCH_OR_f, CBCH_CA_f)

## Zerofill (add zeros to checklists that didn't encounter any CBCH).

# Read in sampling data

df_sampling <- read_sampling("./Data/Raw/EBD_Sampling/sampling_f_relApr-2023.txt", sep = "\t")

df_sampling <- df_sampling %>%
  filter(county %in% CBCH_f$county)

# Zero-fill

CBCH_zf <- auk_zerofill(CBCH_f, df_sampling, collapse=TRUE)

## Apply recommended data-cleaning steps (see 
## https://cornelllabofornithology.github.io/ebird-best-practices/ebird.html):

time_to_decimal <- function(x) {
  x <- hms(x, quiet = TRUE)
  hour(x) + minute(x) / 60 + second(x) / 3600
}

CBCH <- CBCH_zf %>% 
  mutate(
    # convert X to NA
    observation_count = if_else(observation_count == "X", 
                                NA_character_, observation_count),
    observation_count = as.integer(observation_count),
    # effort_distance_km to 0 for non-travelling counts
    effort_distance_km = if_else(protocol_type != "Traveling", 
                                 0, effort_distance_km),
    # convert time to decimal hours since midnight
    time_observations_started = time_to_decimal(time_observations_started),
    # split date into year and day of year
    year = year(observation_date),
    day_of_year = yday(observation_date)
    
  )

## Filter the data to a final analysis package of checklists lasting less than 
## 5 hours, travelling less than 10km, and with 10 or fewer observers, and
## create a "type" variable to facilitate the test-train split for later use
## in model validation

CBCH_zf_f <- CBCH %>% 
  filter(
    # effort filters
    duration_minutes <= 5 * 60,
    effort_distance_km <= 10,
    # 10 or fewer observers
    number_observers <= 10) %>%
  mutate(type = if_else(runif(nrow(.)) <= 0.8, "train", "test"))

# Select desired columns.

CBCH_ebird <- CBCH_zf_f %>% 
  select(checklist_id, observer_id, sampling_event_identifier,
         scientific_name,
         observation_count, species_observed, 
         state_code, county, county_code, locality_id, latitude, 
         longitude, protocol_type, all_species_reported,
         observation_date, year, day_of_year,
         time_observations_started, 
         duration_minutes, effort_distance_km,
         number_observers, type)


## Write

write_csv(CBCH_ebird, "./Data/Cleaned/eBird-Basic-Dataset/Zero-Filled-and-Filtered/CBCH/ebd_CBCH_BCR5_zf.csv")

### Filtering Black-capped Chickadee Data

## Load downloaded eBird Basic Dataset for Washington and sampling database.

BCCH_WA <- auk_ebd("./Data/Raw/eBird-Basic-Dataset/BCCH/ebd_US-WA_bkcchi_relApr-2023.txt", 
                   sep = "\t")

# Filter to BCR 5, complete checklists, and Stationary and Traveling methods.

BCCH_filters_WA <- BCCH_WA %>%
  auk_bcr(bcr = 5) %>%
  auk_protocol(protocol = c("Stationary", "Traveling")) %>%
  auk_complete()

auk_filter(BCCH_filters_WA, 
           file = "./Data/Cleaned/eBird-Basic-Dataset/Filtered/BCCH/WA_BCCH_f_relApr-2023.txt", 
           overwrite = TRUE)

# Load filtered file back into R environment.

BCCH_WA_f <- read_ebd("./Data/Cleaned/eBird-Basic-Dataset/Filtered/BCCH/WA_BCCH_f_relApr-2023.txt",
                      sep = "\t")

## Repeat what was applied to Washington with Oregon, using the filtered sampling
## dataset we created above.

BCCH_OR <- auk_ebd("./Data/Raw/eBird-Basic-Dataset/BCCH/ebd_US-OR_bkcchi_relApr-2023.txt", sep = "\t")

BCCH_filters_OR <- BCCH_OR %>%
  auk_bcr(bcr = 5) %>%
  auk_protocol(protocol = c("Stationary", "Traveling")) %>%
  auk_complete()

auk_filter(BCCH_filters_OR, file = "./Data/Cleaned/eBird-Basic-Dataset/Filtered/BCCH/OR_BCCH_f_relApr-2023.txt",
           overwrite = TRUE)

BCCH_OR_f <- read_ebd("./Data/Cleaned/eBird-Basic-Dataset/Filtered/BCCH/OR_BCCH_f_relApr-2023.txt", 
                      sep = "\t")

## Repeat with British Columbia.

BCCH_BC <- auk_ebd("./Data/Raw/eBird-Basic-Dataset/BCCH/ebd_CA-BC_bkcchi_relApr-2023.txt", sep = "\t")

BCCH_filters_BC <- BCCH_BC %>%
  auk_bcr(bcr = 5) %>%
  auk_protocol(protocol = c("Stationary", "Traveling")) %>%
  auk_complete()

auk_filter(BCCH_filters_BC, file = "./Data/Cleaned/eBird-Basic-Dataset/Filtered/BCCH/BC_BCCH_f_relApr-2023.txt",
           overwrite = TRUE)

BCCH_BC_f <- read_ebd("./Data/Cleaned/eBird-Basic-Dataset/Filtered/BCCH/BC_BCCH_f_relApr-2023.txt", 
                      sep = "\t")

## Repeat with California.

BCCH_CA <- auk_ebd("./Data/Raw/eBird-Basic-Dataset/BCCH/ebd_US-CA_bkcchi_relApr-2023.txt", sep = "\t")

BCCH_filters_CA <- BCCH_CA %>%
  auk_bcr(bcr = c(5, 32)) %>%
  auk_protocol(protocol = c("Stationary", "Traveling")) %>%
  auk_complete()

auk_filter(BCCH_filters_CA, file = "./Data/Cleaned/eBird-Basic-Dataset/Filtered/BCCH/CA_BCCH_f_relApr-2023.txt",
           overwrite = TRUE)

BCCH_CA_f <- read_ebd("./Data/Cleaned/eBird-Basic-Dataset/Filtered/BCCH/CA_BCCH_f_relApr-2023.txt", 
                      sep = "\t")

BCCH_CA_f <- BCCH_CA_f %>%
  filter(county %in% c(unique(CBCH_CA_f$county[CBCH_CA_f$bcr_code == 5]), 
                       "Napa", "Marin", "Solano", "Alameda", "San Mateo", "San Francisco",
                       "Contra Costa", "Santa Clara"))

## Compile into one object

BCCH_f <- rbind(BCCH_BC_f, BCCH_WA_f, BCCH_OR_f, BCCH_CA_f)

## Zerofill (add zeros to checklists that didn't encounter any BCCH).

# Zero-fill

BCCH_zf <- auk_zerofill(BCCH_f, df_sampling, collapse=TRUE)

## Apply recommended data-cleaning steps (see 
## https://cornelllabofornithology.github.io/ebird-best-practices/ebird.html):

BCCH <- BCCH_zf %>% 
  mutate(
    # convert X to NA
    observation_count = if_else(observation_count == "X", 
                                NA_character_, observation_count),
    observation_count = as.integer(observation_count),
    # effort_distance_km to 0 for non-travelling counts
    effort_distance_km = if_else(protocol_type != "Traveling", 
                                 0, effort_distance_km),
    # convert time to decimal hours since midnight
    time_observations_started = time_to_decimal(time_observations_started),
    # split date into year and day of year
    year = year(observation_date),
    day_of_year = yday(observation_date)
    
  )

## Filter the data to a final analysis package of checklists lasting less than 
## 5 hours, travelling less than 10km, and with 10 or fewer observers, and
## create a "type" variable to facilitate the test-train split for later use
## in model validation

BCCH_zf_f <- BCCH %>% 
  filter(
    # effort filters
    duration_minutes <= 5 * 60,
    effort_distance_km <= 10,
    # 10 or fewer observers
    number_observers <= 10) %>%
  mutate(type = if_else(runif(nrow(.)) <= 0.8, "train", "test"))

# Select desired columns.

BCCH_ebird <- BCCH_zf_f %>% 
  select(checklist_id, observer_id, sampling_event_identifier,
         scientific_name,
         observation_count, species_observed, 
         state_code, county, county_code, locality_id, latitude, 
         longitude, protocol_type, all_species_reported,
         observation_date, year, day_of_year,
         time_observations_started, 
         duration_minutes, effort_distance_km,
         number_observers, type)


## Write

write_csv(BCCH_ebird, "./Data/Cleaned/eBird-Basic-Dataset/Zero-Filled-and-Filtered/BCCH/ebd_BCCH_BCR5_zf.csv")
