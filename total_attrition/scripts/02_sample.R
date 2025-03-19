
# Script for sampling coefficients and variance-covariance matrices from a
# multivariate normal distribution

source("config.R")
library(tidyverse)

# Import prepared data
lst_prep <- read_rds("total_attrition/processed_data/prepared_data.rds")

# Check Kaplan-Meier curves
plot_cind <- unnest(lst_prep$cind, cols = "t2e") %>% 
  ggplot(aes(x = time, y = 1 - estimate)) +
  geom_step(linewidth = 1) +
  facet_wrap(~ctgov, scales = "free") +
  labs(x = "Time in trial (days", y = "Cumulative incidence of attrition")