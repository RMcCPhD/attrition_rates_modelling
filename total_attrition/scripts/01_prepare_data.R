
# Script for preparing data exported from vivli for public download and analysis
# Public versions are saved as csv for easier access

source("config.R")
library(tidyverse)

### Empirical cumulative distribution ##########################################

# Import exported version
imp_cind <- read_csv("total_attrition/data/t2e.csv")

# Separate variables
# Data for 92 trials extracted but two were removed when finalising the
# analysis set:
# NCT01131676 - Discrepancies noted between the extracted counts and ctgov
# NCT00274573 - Trial was stopped by the sponsor, not representative of attr
cind_sep <- imp_cind %>% 
  filter(!(ctgov %in% c("NCT01131676", "NCT00274573"))) %>% 
  group_by(ctgov) %>% 
  separate_rows(dat, sep = " ; ") %>% 
  separate(dat, c("time", "estimate"), sep = ", ") %>% 
  ungroup() %>% 
  mutate(
    across(c(time, estimate), ~ gsub("t|est|\\=", "", .)),
    across(c(time, estimate), ~ as.numeric(.))
  )

glimpse(cind_sep)

# Save as csv for public access
# write_csv(cind_sep, "usable_data/total_attrition/cind_pblc.csv")

# Save as nested data for analysis
cind_nest <- cind_sep %>% 
  group_by(ctgov) %>% 
  nest(.key = "t2e") %>% 
  ungroup()

### Best parametric models #####################################################

# Import exported versions
imp_params <- read_csv("total_attrition/data/parameters.csv")
imp_vcov <- read_csv("total_attrition/data/vcov.csv")
imp_fit <- read_csv("total_attrition/data/fit.csv")

# Get the best fitting models for each trial
# Convergence issues with generalised gamma so omitted from selections
mdl_best <- imp_fit %>% 
  filter(dist != "gengamma" & !(ctgov %in% c("NCT01131676", "NCT00274573"))) %>% 
  group_by(ctgov) %>% 
  slice_min(aic, n = 1) %>% 
  slice_head(n = 1) %>% 
  ungroup()

# Join parameters for the best-fitting models
# Nest parameters
mdl_params_fit <- mdl_best %>% 
  left_join(imp_params) %>% 
  group_by(ctgov, loglik, aic) %>% 
  nest(.key = "params") %>% 
  ungroup()

# Get variance-covariance matrices for best-fitting models
# Nest
mdl_vcov <- imp_vcov %>% 
  inner_join(mdl_best %>% select(ctgov, dist) %>% distinct()) %>% 
  group_by(ctgov) %>% 
  nest(.key = "vcov") %>% 
  ungroup()

# Create a public version saved as csv files for readability
mdl_params_pblc <- unnest(mdl_params_fit, cols = "params")

mdl_vcov_pblc <- unnest(mdl_vcov, cols = "vcov") %>% 
  rename(param1 = parameter1, param2 = parameter2) %>% 
  mutate(param1 = if_else(param1 == "1", "rate", param1))

# write_csv(mdl_params_pblc, "usable_data/total_attrition/params_pblc.csv")
# write_csv(mdl_vcov_pblc, "usable_data/total_attrition/vcov_pblc.csv")

# Convert nested data in parameters and vcov to a format usable for sampling
# from a multivariate normal distribution
mdl_fmt_params <- mdl_params_fit %>% 
  mutate(
    params = map(params, function(df) {
      df %>%
        select(-se) %>%
        group_by(dist) %>%
        pivot_wider(names_from = parameter, values_from = est) %>%
        nest(params = -dist) %>%
        arrange(dist) %>%
        mutate(params = map(params, ~ as.matrix(select(.x, where(~ all(!is.na(.)))))))
    })
  )

mdl_fmt_vcov <- mdl_vcov %>% 
  mutate(
    vcov = map(vcov, function(df) {
      df %>%
        group_by(dist) %>%
        pivot_wider(names_from = parameter2, values_from = r) %>%
        nest(vcov = -dist) %>%
        ungroup() %>%
        mutate(
          vcov = map(vcov, ~ {
            mat <- .x %>%
              select(where(~ all(!is.na(.)))) %>%
              select(-parameter1) %>%
              as.matrix()
            
            rownames(mat) <- colnames(mat)
            mat
          })
        )
    })
  ) %>% 
  arrange(ctgov)

# Save as a list of prepared data for analysis
lst_dat <- list(cind = cind_nest, params = mdl_fmt_params, vcov = mdl_fmt_vcov)
# saveRDS(lst_dat, "total_attrition/processed_data/prepared_data.rds")
