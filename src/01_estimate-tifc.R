#-----------------------------------------------------------------------------------------------------------------------

# Title:       Estimate time in front of camera for paired heights analysis
# Description: At each of the 10 paired sites with two cameras at different heights (1m vs 0.5m), estimate the total
#              time spent in front of the camera for each mammal species.

# Author:      Marcus A Becker
# Date:        July 28, 2022

#-----------------------------------------------------------------------------------------------------------------------

# Root directory (Shared Google Drive)
root <- "G:/Shared drives/ABMI Camera Mammals/"

# Attach packages
library(wildRtrax) # Just obtaining data
library(keyring) # Storing credentials

library(dplyr)
library(stringr)
library(tidyr)
library(purrr)
library(lubridate)
library(readr)

# Native species tags in WildTrax
load(paste0(root, "data/lookup/wt_native_sp.RData"))

# Probabilistic gaps probability of leaving predictions
df_leave_prob_pred <- read_csv(paste0(root, "data/processed/probabilistic-gaps/gap-leave-prob_predictions_2021-10-05.csv"))
# Gap groups
df_gap_groups <- read_csv(paste0(root, "data/lookup/species-gap-groups.csv"))
# Average time between images per species
df_tbp <- read_csv(paste0(root, "data/processed/time-btwn-images/abmi-cmu_all-years_tbp_2021-06-25.csv"))

# Seasonal start/end dates (julian day):
summer.start.j <- 106 # April 16
summer.end.j <- 288 # October 15

#-----------------------------------------------------------------------------------------------------------------------

# Load data

# Set up credentials for WildTrax
Sys.setenv(WT_USERNAME = key_get("WT_USERNAME", keyring = "wildtrax"),
           WT_PASSWORD = key_get("WT_PASSWORD", keyring = "wildtrax"))

# Authenticate
wt_auth()

# Obtain database project IDs
proj_ids <- wt_get_download_summary(sensor_id = "CAM") |>
  filter(str_detect(project, "Ecosystem Health 2022|Height")) |>
  select(project_id) |>
  pull() |>
  unlist()

# Deployments of interest (to filter out from EH 2022)
locations <- c("759-NE", "759-SW", "760-NE", "760-SW", "792-NE",
               "792-SW", "793-NE", "793-SW", "794-NE", "794-SW")

# In Ecosystem Health 2022 (1m cameras):
# - 793-NE failed August 6, 2021
# - 793-SW failed May 18, 2022

# In Heights (0.5m cameras):
# - 792-SW failed July 1, 2022

# Download data
data <- map_df(.x = proj_ids,
               .f = ~ wt_download_report(
                 project_id = .x,
                 sensor_id = "CAM",
                 cols_def = FALSE,
                 weather_cols = FALSE
               )) |>
  mutate(date_detected = ymd_hms(date_detected)) |>
  filter(location %in% locations,
         # Truncate dates so that only data operating during a common period is used (only 3/10 locations impacted)
         !(date_detected > as.Date("2021-08-07 00:00:00") & location == "793-NE"),
         !(date_detected > as.Date("2022-05-19 00:00:00") & location == "793-SW"),
         !(date_detected > as.Date("2022-07-02 00:00:00") & location == "792-SW")) |>
         select(project, location, date_detected, field_of_view, common_name, number_individuals)

#-----------------------------------------------------------------------------------------------------------------------

# Assess where NONE images split detections
df_gap_nones <- data %>%
  select(project, location, date_detected, common_name) %>%
  arrange(project, location, date_detected) %>%
  # Create gap class column
  mutate(common_name_next = lead(common_name),
         gap_class = ifelse(common_name != "NONE" & common_name_next == "NONE", "N", NA)) %>%
  filter(gap_class == "N") %>%
  select(-c(common_name_next))

# Identify independent detections ('series')
df_series <- data |>
  filter(common_name %in% native_sp,
         field_of_view == "WITHIN") |>
  # Identify where NONEs occurred
  left_join(df_gap_nones, by = c("location", "project", "date_detected", "common_name")) |>
  # Order observations
  arrange(project, location, date_detected, common_name) |>
  # Identify series
  mutate(series_num = 0,
         # Lagged time stamp
         date_detected_previous = lag(date_detected),
         # Lead time stamp
         date_detected_next = lead(date_detected),
         # Calculate difference in time between ordered images
         diff_time_previous = as.numeric(date_detected - date_detected_previous),
         diff_time_next = as.numeric(abs(date_detected - date_detected_next)),
         # Lagged species
         common_name_previous = lag(common_name),
         # Was is a different species?
         diff_sp = ifelse(common_name != common_name_previous, TRUE, FALSE),
         # Lagged deployment
         location_previous = lag(location),
         # Was is a different deployment?
         diff_location = ifelse(location != location_previous, TRUE, FALSE),
         # Flag gaps that will need checking
         gap_check = ifelse(diff_location == FALSE & diff_sp == FALSE & (diff_time_previous <= 120 & diff_time_previous >= 20), 1, 0),
         # Lagged gap class
         gap_class_previous = replace_na(lag(gap_class), ""),
         # Identify new series, based on being different deployment, species, greater than 120 seconds, and approp gaps
         diff_series = ifelse(diff_location == TRUE | diff_sp == TRUE | diff_time_previous > 120 | (gap_class_previous == "L" | gap_class_previous == "N"), 1, 0),
         # Number series
         series_num = c(0, cumsum(diff_series[-1])),
         # Flag gaps that require probabilistic time assignment
         gap_prob = replace_na(ifelse(gap_check == 1 & (gap_class_previous == "" | gap_class_previous == "U"), 1, 0), 0)) |>
  group_by(series_num) |>
  mutate(diff_time_previous = ifelse(row_number() == 1, 0, diff_time_previous),
         diff_time_next = ifelse(row_number() == n(), 0, diff_time_next)) |>
  ungroup() |>
  # Join gap group lookup table
  left_join(df_gap_groups, by = "common_name") |>
  # Join gap leaving predictions
  left_join(df_leave_prob_pred, by = c("gap_group", "diff_time_previous" = "diff_time")) |>
  # Adjust time difference between ordered images that require probabilistic time assignment
  mutate(pred = replace_na(pred, 1),
         diff_time_previous_adj = ifelse(gap_prob == 1, diff_time_previous * (1 - pred), diff_time_previous),
         diff_time_next_adj = ifelse(lead(gap_prob == 1), diff_time_next * (1 - lead(pred)), diff_time_next))

# Vector of all project-locations
dep <- df_series |>
  unite(col = "location_project", location, project, sep = "_") |>
  select(location_project) |>
  distinct() |>
  pull()

# Total time for each detection/series
df_tts <- df_series |>
  left_join(df_tbp, by = "common_name") |>
  filter(!number_individuals == "VNA") |>
  mutate(number_individuals = as.numeric(number_individuals)) |>
  group_by(series_num) |>
  mutate(# Check whether the image was first or last in a series
    bookend = ifelse(row_number() == 1 | row_number() == n(), 1, 0),
    # Calculate time for each individual image
    image_time = ifelse(bookend == 1,
                        ((diff_time_previous_adj + diff_time_next_adj) / 2) + (tbp / 2),
                        (diff_time_previous_adj + diff_time_next_adj) / 2),
    # Multiply image time by the number of animals present
    image_time_ni = image_time * number_individuals) |>
  # Group by common name as well to add it as a variable to output
  group_by(common_name, location, project, .add = TRUE) |>
  # Calculate total time and number of images for each series
  summarise(n_images = n(),
            series_total_time = sum(image_time_ni)) |>
  ungroup() |>
  # Just adjust a single WTD series
  mutate(series_total_time = if_else(is.na(series_total_time), 22, series_total_time)) |>
  select(series_num, common_name, n_images, series_total_time)

# Calculate total time in front of the camera, by deployment, project, and species
df_tt <- df_series |>
  group_by(series_num) |>
  arrange(date_detected, .by_group = TRUE) |>
  filter(row_number() == 1) |>
  left_join(df_tts, by = c("series_num", "common_name")) |>
  select(project, location, date_detected, common_name, series_num, series_total_time) |>
  ungroup() |>
  mutate(julian = as.numeric(format(date_detected, "%j")),
         season = ifelse(julian >= summer.start.j & julian <= summer.end.j, "summer", "winter")) |>
  unite(location, project, col = "location_project", sep = "_", remove = TRUE) |>
  mutate_at(c("location_project", "common_name", "season"), factor) |>
  # Note: didn't include season here - not enough sample size at this point (I think ...). TBD.
  group_by(location_project, common_name, .drop = FALSE) |>
  summarise(total_duration = sum(series_total_time)) |>
  ungroup() |>
  mutate_if(is.factor, as.character)

#-----------------------------------------------------------------------------------------------------------------------

# Save results
write_csv(df_tt, "./data/processed/heights-experiment_tifc_long.csv")

#-----------------------------------------------------------------------------------------------------------------------
