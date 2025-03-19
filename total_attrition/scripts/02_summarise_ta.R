
# Load packages, functions, etc.
source("Scripts/00b_packages_functions_total.R")

# Import data
import_covariance <- read_csv("../Analysis/Vivli/Models_ta/covariance.csv")
import_fit <- read_csv("../Analysis/Vivli/Models_ta/fit.csv")
import_parameters <- read_csv("../Analysis/Vivli/Models_ta/parameters.csv")
import_t2e <- read_csv("../Analysis/Vivli/Models_ta/t2e.csv")

data <- read_rds("../Analysis/Data/ta/data.rds")
sample_parameters <- read_rds("../Analysis/Data/ta/sample_parameters.rds")
sample_estimates <- read_rds("../Analysis/Data/ta/sample_estimates.rds")

import_flow <- read_csv("../Analysis/Vivli/Summary_statistics/participant_flow_formatted.csv")

# Create vector to order conditions
order_conditions <- c(
  "Type 2 diabetes", "COPD", "Hypertension", "Benign prostatic hyperplasia", 
  "Parkinson's disease", "Pulmonary fibrosis", "Asthma", "Osteoarthritis", 
  "Psoriasis", "Restless leg syndrome", "Osteoporosis"
)

import_summary <- import_flow %>% 
  mutate(
    condition = factor(condition, levels = order_conditions),
    across(
      completers:dropouts,
      ~ gsub("\\s*\\([^)]*\\)", "", .)
    ),
    across(
      completers:dropouts,
      as.numeric
    )
  ) %>% 
  arrange(condition)

# Check frequency of best fitting models
tbl_best_models <- unnest(data$best_models, c("nested_params")) %>%
  inner_join(import_summary %>% dplyr::select(ctgov)) %>% 
  dplyr::select(ctgov, dist) %>% 
  distinct() %>% 
  group_by(dist) %>% 
  reframe(n = n()) %>% 
  mutate(
    pcnt = round((n/sum(n))*100, 1),
    values = paste0(n , " (", pcnt, ")"),
    total = sum(n),
    total_pcnt = round(total/sum(total) * 100, 1),
    total_value = paste0(total , " (", total_pcnt, ")")
  ) %>% 
  pivot_wider(
    names_from = dist,
    values_from = values
  ) %>% 
  dplyr::select(
    contains(c("exp")),
    # contains(c("gengamma")),
    contains(c("gompertz")),
    contains(c("lnorm")),
    contains(c("llogis")),
    contains(c("weibull")),
    contains(c("total_value"))
  ) %>% 
  fill(everything(), .direction = "updown") %>% 
  ungroup() %>% 
  distinct()

write_csv(
  tbl_best_models,
  "../Analysis/Output/ta_tbls/best_models.csv"
)

## Comparing Gompertz and log-normal -------------------------------------------

# Get table comparing AIC scores for best models
compare_aic <- unnest(data$best_models, c("nested_params")) %>%
  inner_join(import_summary %>% dplyr::select(ctgov)) %>% 
  dplyr::select(ctgov, dist) |> 
  filter(dist %in% c("lnorm", "gompertz")) |> 
  dplyr::select(ctgov) |> 
  inner_join(
    data$parameters_fit |> 
      dplyr::select(ctgov, dist, aic) |> 
      distinct()
  ) |> 
  filter(dist %in% c("lnorm", "gompertz")) |> 
  pivot_wider(
    names_from = dist,
    values_from = aic
  ) |> 
  mutate(
    aic_diff = abs(gompertz - lnorm)
  )

# Summarise comparison of AIC between gompertz and log-normal
# Median: 3
summary(compare_aic$aic_diff)

# Get trials within quantiles (0-7)
quantile_trials <- compare_aic |> filter(aic_diff < 7)

# Get model parameters for these trials
parameters_trials <- quantile_trials |> 
  pivot_longer(
    gompertz:lnorm,
    names_to = "dist",
    values_to = "aic"
  ) |> 
  dplyr::select(ctgov, dist) |> 
  inner_join(
    data$parameters_fit |> 
      dplyr::select(-loglik, -aic, -se)
  )

# Get KM estimates for these trials
km_trials <- quantile_trials |> 
  dplyr::select(ctgov) |> 
  inner_join(import_t2e) |> 
  group_by(ctgov) %>% 
    separate_rows(
      dat,
      sep = " ; "
    ) %>% 
    separate(
      dat,
      into = c("time", "estimate"),
      sep = ", "
    ) %>% 
    ungroup() %>% 
    mutate(
      across(
        c(time, estimate),
        ~ gsub("t|est|\\=", "", .)
      ),
      across(
        c(time, estimate),
        ~ as.numeric(.)
      ),
      estimate = 1 - estimate
    )

# Sequence times for each trial
sequence_km_times <- km_trials |> 
  dplyr::select(ctgov, time) |> 
  group_by(ctgov) |> 
  summarise(
    min_time = min(time),
    max_time = max(time)
  ) |> 
  mutate(
    time = map2(min_time, max_time, seq)
  ) |> 
  unnest(cols = time) |> 
  dplyr::select(ctgov, time) |> 
  left_join(km_trials)

# Estimate cdf at each time point per trial (Gompertz)
estimate_cdf_gompertz <- sequence_km_times |> 
  left_join(
    parameters_trials |> 
      filter(dist == "gompertz") |> 
      dplyr::select(ctgov, parameter, est) |> 
      pivot_wider(
        names_from = "parameter",
        values_from = "est"
      )
  ) |> 
  group_by(ctgov, time) |> 
  mutate(
    mdl_cdf = pgompertz(q = time, shape = shape, rate = exp(rate))
  ) |> 
  ungroup() |> 
  dplyr::select(-rate, -shape) |> 
  rename(
    km = estimate,
    gompertz = mdl_cdf
  ) |> 
  pivot_longer(
    km:gompertz,
    names_to = "type",
    values_to = "estimate"
  ) |> 
  na.omit()

# Estimate cdf at each time point per trial (log-normal)
estimate_cdf_lnorm <- sequence_km_times |> 
  left_join(
    parameters_trials |> 
      filter(dist == "lnorm") |> 
      dplyr::select(ctgov, parameter, est) |> 
      pivot_wider(
        names_from = "parameter",
        values_from = "est"
      )
  ) |> 
  group_by(ctgov, time) |> 
  mutate(
    mdl_cdf = plnorm(q = time, meanlog = meanlog, sdlog = exp(sdlog))
  ) |> 
  ungroup() |> 
  dplyr::select(-meanlog, -sdlog) |> 
  rename(
    km = estimate,
    lnorm = mdl_cdf
  ) |> 
  pivot_longer(
    km:lnorm,
    names_to = "type",
    values_to = "estimate"
  ) |> 
  na.omit()

# Combine both into one dataframe
estimate_cdf <- estimate_cdf_gompertz |> 
  full_join(estimate_cdf_lnorm) |> 
  arrange(ctgov, time) |> 
  pivot_wider(
    names_from = "type",
    values_from = "estimate"
  )

# Plot estimates against KM CDF
plot_cdf <- estimate_cdf |> 
  # mutate(
  #   new_id = dense_rank(ctgov),
  #   new_id = if_else(new_id < 10L, paste0("0", new_id), as.character(new_id)),
  #   new_id = paste0("Trial ", new_id),
  #   new_id = factor(new_id)
  # ) |>
  filter(!is.na(km)) |> 
  ggplot(
    aes(
      x = time,
      y = km
    )
  ) +
  geom_step(
    linewidth = 0.5,
    colour = "black"
  ) +
  geom_line(
    aes(
      x = time,
      y = gompertz
    ),
    colour = "orange",
    linewidth = 0.5
  ) +
  geom_line(
    aes(
      x = time,
      y = lnorm
    ),
    colour = "purple",
    linewidth = 0.5
  ) +
  scale_x_continuous(n.breaks = 4) +
  scale_y_continuous(n.breaks = 4) +
  labs(
    x = "Time in days",
    y = "Cumulative incidence of attrition (%)"
  ) +
  facet_wrap(~ctgov, scales = "free") +
  theme_bw() +
  theme(
    text = element_text(size = 9),
  )

ggsave(
  "../Analysis/Output/plots/gompertz_lnorm_compare.png",
  plot_cdf,
  dpi = 300,
  width = 10,
  height = 6
)

## Summarise hazard rates ------------------------------------------------------

# Import tally
import_tally_hazards <- read_csv("../Analysis/Output/ta_tbls/tally_hazards_total.csv") %>% 
  inner_join(import_summary %>% dplyr::select(ctgov))
  
# High point
tbl_hazard_highest <- import_tally_hazards %>% 
  group_by(dist, peak) %>% 
  reframe(n = n()) %>% 
  group_by(dist) %>% 
  mutate(
    sum_n = sum(n),
    pcnt = round((n/sum_n) * 100, 1),
    value = paste0(n, " (", pcnt, ")")
  ) %>% 
  ungroup() %>% 
  pivot_wider(
    names_from = "peak",
    values_from = "value"
  ) %>% 
  write_csv(
    "../Analysis/Output/ta_tbls/hazards_peak.csv"
  )

# Model and direction
tbl_hazard_model_direction <- import_tally_hazards %>% 
  group_by(dist, direction) %>% 
  reframe(n = n()) %>% 
  group_by(dist) %>% 
  mutate(
    sum_n = sum(n),
    pcnt = round((n/sum_n) * 100, 1),
    value = paste0(n, " (", pcnt, ")")
  ) %>% 
  ungroup() %>% 
  pivot_wider(
    names_from = "direction",
    values_from = "value"
  ) %>% 
  write_csv(
    "../Analysis/Output/ta_tbls/hazards_direction.csv"
  )

# Model and shape
tbl_hazard_model_shape <- import_tally_hazards %>% 
  group_by(dist, shape) %>% 
  reframe(n = n()) %>% 
  group_by(dist) %>% 
  mutate(
    sum_n = sum(n),
    pcnt = round((n/sum_n) * 100, 1),
    value = paste0(n, " (", pcnt, ")")
  ) %>% 
  ungroup() %>% 
  pivot_wider(
    names_from = "shape",
    values_from = "value"
  ) %>% 
  write_csv(
    "../Analysis/Output/ta_tbls/hazards_shape.csv"
  )









