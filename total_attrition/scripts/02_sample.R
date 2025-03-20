
# Script for sampling coefficients and variance-covariance matrices from a
# multivariate normal distribution and producing datasets for plotting
# Outputs are not saved in github repo as file sizes are quite large, but
# the script can be run using uploaded data to get the outputs

source("config.R")
source("total_attrition/scripts/00_functions.R")

library(tidyverse)
library(MASS)
library(flexsurv)

# Import prepared data
lst_prep <- read_rds("total_attrition/processed_data/prepared_data.rds")

# Check empirical cumulative incidence curves
plot_cind <- unnest(lst_prep$cind, cols = "cind") %>% 
  ggplot(aes(x = time, y = 1 - estimate)) +
  geom_step(linewidth = 1) +
  facet_wrap(~ctgov, scales = "free") +
  labs(x = "Time in trial (days)", y = "Cumulative incidence of attrition")

plot_cind

# Sample parameters from a multivariate distribution for each trial model
sample_parameters <- get_cind()

# Estimate cumulative incidence using sampled parameters
# Takes about an hour to run
cind_sampled_estimates <- get_estimates()

# Save this list to prevent having to re-run (quite a large file so have not
# added this to the github repo, but the script can be run using the provided
# data to create this)
# saveRDS(cind_sampled_estimates, "total_attrition/processed_data/cind_sample_estimates.rds")

# Merge cumulative incidence datasets
cind_bind <- lapply(seq_along(cind_sampled_estimates), function(i) {
  name <- names(cind_sampled_estimates[i])
  lapply(seq_along(cind_sampled_estimates[[i]]), function(j) {
    name2 <- names(cind_sampled_estimates[[i]][j])
    bind_rows(cind_sampled_estimates[[i]][[j]]) %>% 
      mutate(
        ctgov = name,
        dist = name2
      )
  })
})

cind_bind <- bind_rows(cind_bind) %>% distinct()
# saveRDS(cind_bind, "total_attrition/processed_data/cind_sampled_df.rds")

# Estimate hazard rates over time using model parameters
haz_bind <- get_hazards()

# Save (locally, can run the script as said above to get these objects)
# saveRDS(haz_bind, "total_attrition/processed_data/haz_sampled_df.rds")

# Unload MASS after use to avoid dplyr conflicts
detach("package:MASS", unload = TRUE)