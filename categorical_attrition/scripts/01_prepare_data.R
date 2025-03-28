
# Script for preparing data exported from vivli for public download and analysis
# Public versions are saved as csv for easier access

source("config.R")
library(tidyverse)

### Empirical cumulative distribution ##########################################

# Import exported version
imp_cind <- read_csv("categorical_attrition/data/t2e.csv")

# Note of amendments
# See total_attrition/01_prepare_data.R
cind_sep <- imp_cind %>% 
  filter(!(ctgov %in% c("NCT01131676", "NCT00274573"))) %>% 
  group_by(ctgov) %>% 
  separate_rows(dat, sep = " ; ") %>% 
  separate(dat, c("time", "estimate"), sep = ", ") %>% 
  ungroup() %>% 
  mutate(
    across(c(time, estimate), ~ gsub("t|est|\\=", "", .)),
    across(c(time, estimate), ~ as.numeric(.))
  ) %>% 
  group_by(ctgov) %>% 
  mutate(
    estimate = case_when(
      ctgov == "NCT01694771" 
      & time > 84
      & !is.na(estimate) ~ first(estimate[time == 84]),
      TRUE ~ estimate
    )
  ) %>% 
  ungroup() %>% 
  filter(
    !(ctgov == "NCT01694771" & time >= 84),
    !(ctgov == "NCT00384930" & time >= 84),
    !(ctgov == "NCT01306214" & time >= 365),
    !(ctgov == "NCT01422876" & time >= 365)
  )

glimpse(cind_sep)

# Test plot
# cind_sep %>% 
#   ggplot(aes(x = time, y = 1 - estimate, colour = cause)) +
#   geom_line(linewidth = 1) +
#   facet_wrap(~ctgov, scales = "free")

# Save as csv for public access
# write_csv(cind_sep, "usable_data/categorical_attrition/cind_pblc.csv")

# Save as nested data for analysis
cind_nest <- cind_sep %>% 
  group_by(ctgov, cause) %>% 
  nest(.key = "cind") %>% 
  ungroup()

### Best parametric models #####################################################

# Import exported versions
imp_params <- read_csv("categorical_attrition/data/parameters.csv")
imp_vcov <- read_csv("categorical_attrition/data/vcov.csv")
imp_fit <- read_csv("categorical_attrition/data/fit.csv")

# Get the best fitting models for each trial
# Convergence issues with generalised gamma so omitted from selections
mdl_best <- imp_fit %>% 
  filter(dist != "gengamma" & !(ctgov %in% c("NCT01131676", "NCT00274573"))) %>% 
  group_by(ctgov, cause) %>% 
  slice_min(aic, n = 1) %>% 
  ungroup()

# Where multiple models had the lowest AIC (22 trials), select the one that
# best fits empirical cumulative incidence
mdl_best2 <- mdl_best %>% 
  left_join(
    read_csv("vivli/scratch_data/n22_log.csv") %>% distinct()
  ) %>% 
  group_by(ctgov, cause) %>% 
  mutate(
    take = case_when(
      n() == 1 ~ 1,
      n() != 1 & !is.na(best_fit) & best_fit == dist ~ 1,
      TRUE ~ 0
    )
  ) %>% 
  ungroup() %>% 
  filter(take == 1) %>% 
  select(-take, -best_fit)

# Join parameters for the best-fitting models
# Nest parameters
mdl_params_fit <- mdl_best2 %>% 
  left_join(imp_params) %>% 
  arrange(ctgov, cause) %>% 
  group_by(ctgov, cause, loglik, aic) %>% 
  nest(.key = "params") %>% 
  ungroup()

# Get variance-covariance matrices for best-fitting models
# Nest
mdl_vcov <- imp_vcov %>% 
  inner_join(mdl_best2 %>% select(ctgov, cause, dist) %>% distinct()) %>% 
  arrange(ctgov, cause) %>% 
  group_by(ctgov, cause) %>% 
  nest(.key = "vcov") %>% 
  ungroup()

# Create a public version saved as csv files for readability
mdl_params_pblc <- unnest(mdl_params_fit, cols = "params")

mdl_vcov_pblc <- unnest(mdl_vcov, cols = "vcov") %>% 
  rename(param1 = parameter1, param2 = parameter2) %>% 
  mutate(param1 = if_else(param1 == "1", "rate", param1))

# write_csv(mdl_params_pblc, "usable_data/categorical_attrition/params_pblc.csv")
# write_csv(mdl_vcov_pblc, "usable_data/categorical_attrition/vcov_pblc.csv")

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
# saveRDS(lst_dat, "categorical_attrition/processed_data/prepared_data.rds")
