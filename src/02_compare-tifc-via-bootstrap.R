#-----------------------------------------------------------------------------------------------------------------------

# Title:       Compare TIFC for each species via bootstrapping
# Description: Compare the calculated time spent in front of the camera between the two height conditions via
#              bootstrapping.

# Author:      Marcus A Becker
# Date:        August 17, 2022

#-----------------------------------------------------------------------------------------------------------------------

# Attach packages
library(readr)
library(tidyr)
library(dplyr)
library(stringr)

# Read in previously processed TIFC data from script 01
df_tt <- read_csv("./data/processed/heights-experiment_tifc_long.csv")

# Only interested in a subset of species with sufficient sample size
sp_of_interest <- c("Black Bear", "White-tailed Deer", "Moose", "Coyote", "Cougar", "Marten",
                    "Snowshoe Hare", "Fisher", "Canada Lynx", "Gray Wolf", "Red Squirrel")

#-----------------------------------------------------------------------------------------------------------------------

# Prepare data for bootstrapping
boot <- df_tt %>%
  filter(common_name %in% sp_of_interest) %>%
  separate(location_project, into = c("location", "project"), sep = "_") %>%
  mutate(height = ifelse(str_detect(project, "Ecosystem"), "full_height", "half_height")) %>%
  pivot_wider(id_cols = c(location, common_name),
              names_from = "height",
              values_from = c("total_duration")) %>%
  mutate(loc_sp = paste0(location, "_", common_name),
         sp = paste0(common_name)) %>%
  select(loc_sp, 3:6, location, sp) %>%
  # Remove deployment-species-season combos without any duration/images at either camera
  filter(full_height > 0 | half_height > 0)

#-----------------------------------------------------------------------------------------------------------------------

# Bootstrap

sp.list <- sort(unique(boot$sp))

niter <- 1000
bs1 <- array(0, c(length(sp.list), 2, niter))

dimnames(bs1)[[1]] <- sp.list

bs.sum<-array(0, c(length(sp.list), 1, 3))
dimnames(bs.sum) <- list(sp.list, c("duration_full_as_pct_of_half"), c("median","q5","q95"))

# Set seed for reproducibility:
set.seed(12345)

for (sp in 1:length(sp.list)) {

  print(paste(sp,length(sp.list),date()))
  boot.sp<-boot[boot$sp==sp.list[sp],]
  loc.list <- sort(unique(boot.sp$location))

  for (iter in 1:niter) {

  # Sample from the locations, with replacement.
  s <- sample(1:length(loc.list), length(loc.list), replace=TRUE)

  boot1<-NULL

    for (j in 1:length(s))

    boot1<-rbind(boot1, boot.sp[boot.sp$location == loc.list[s[j]],])

    x <- mean(boot1$full_height)
    bs1[sp.list[sp], 1, iter] <- as.numeric(x)

    x <- mean(boot1$half_height)
    bs1[sp.list[sp], 2, iter] <- as.numeric(x)

  }

  # Calculate ratio between the values of both cameras
  ratio1<-bs1[sp,2,]/bs1[sp,1,]*100
  ratio1<-ratio1[!is.na(ratio1)]
  # Obtain median, 5% and 95% quantiles for CIs
  bs.sum[sp,1,]<-quantile(ratio1,c(0.5,0.05,0.95))

}

# Collect values from matrices into table
table <- data.frame(
  sp = sort(unique(boot$sp)),
  npairs = as.numeric(by(boot$location, boot$sp, length)),
  total_duration_full = as.numeric(by(boot$full_height, boot$sp, sum)),
  total_duration_half = as.numeric(by(boot$half_height, boot$sp, sum)),
  half_dur_as_pct_of_full_median = bs.sum[,1,1],
  half_dur_as_pct_of_full_lci = bs.sum[,1,2],
  half_dur_as_pct_of_full_uci = bs.sum[,1,3]) |>
  # Round
  mutate(across(.cols = total_duration_full:last_col(), .f = ~ round(., digits = 2)))

#-----------------------------------------------------------------------------------------------------------------------

# Save results

write_csv(table, "./data/processed/heights-experiment_tifc-bootstrap-comp.csv")

#-----------------------------------------------------------------------------------------------------------------------
