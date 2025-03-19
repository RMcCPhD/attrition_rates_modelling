
# Load packages, functions, etc.
source("Scripts/00a_packages_functions_coa.R")

# Import data
import_covariance <- read_csv("../Analysis/Vivli/Models_coa/covariance.csv")
import_fit <- read_csv("../Analysis/Vivli/Models_coa/fit.csv")
import_parameters <- read_csv("../Analysis/Vivli/Models_coa/parameters.csv")
import_t2e <- read_csv("../Analysis/Vivli/Models_coa/t2e.csv")

data <- read_rds("../Analysis/Data/coa/data.rds")
sample_parameters <- read_rds("../Analysis/Data/coa/sample_parameters.rds")
sample_estimates <- read_rds("../Analysis/Data/coa/sample_estimates.rds")

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
      completers:`Voluntary Withdrawal`,
      ~ gsub("\\s*\\([^)]*\\)", "", .)
    ),
    across(
      completers:`Voluntary Withdrawal`,
      as.numeric
    )
  ) %>% 
  arrange(condition)

## Heatmap of model fit --------------------------------------------------------

# Check frequency of best fitting models
tbl_best_models <- unnest(data$best_models, c("nested_params")) %>% 
  inner_join(import_summary %>% dplyr::select(ctgov)) %>% 
  dplyr::select(ctgov, cause, dist) %>% 
  distinct() %>% 
  group_by(cause, dist) %>% 
  reframe(n = n()) %>% 
  group_by(cause) %>% 
  mutate(
    pcnt = round((n/sum(n))*100, 1),
    values = paste0(n , " (", pcnt, ")"),
    total = sum(n),
    total_pcnt = round(total/263 * 100, 1),
    total_value = paste0(total , " (", total_pcnt, ")")
  ) %>% 
  pivot_wider(
    names_from = dist,
    values_from = values
  ) %>% 
  dplyr::select(
    cause, 
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
  "../Analysis/Output/coa_tbls/best_models_coa.csv"
)

# Dataframe for checking condition placement
df_check_placement <- import_summary %>% 
  dplyr::select(ctgov, condition) %>% 
  distinct() %>% 
  mutate(
    condition = factor(condition, levels = order_conditions)
  ) %>% 
  arrange(condition)

# Get heatmap data for frequency of best performing models across causes
heatmap_data <- unnest(data$best_models, c("nested_params")) %>% 
  inner_join(import_summary %>% dplyr::select(ctgov)) %>% 
  dplyr::select(ctgov, cause, dist) %>% 
  distinct() %>% 
  mutate(
    dist = case_when(
      dist == "exp" ~ "Exponential",
      # dist == "gengamma" ~ "Generalised gamma", # No gengamma selected
      dist == "gompertz" ~ "Gompertz",
      dist == "llogis" ~ "Log-logistic",
      dist == "lnorm" ~ "Log-normal",
      dist == "weibull" ~ "Weibull"
    )
  ) %>% 
  complete(ctgov, cause) %>% 
  left_join(
    import_summary %>% 
      dplyr::select(ctgov, condition) %>% 
      distinct()
  ) %>% 
  mutate(
    condition = factor(condition, levels = order_conditions)
  ) %>% 
  arrange(condition)

# Check overall frequency of best performing models
table(heatmap_data$dist)

# Variable for ordering trial IDs on heatmap
var_order_trials <- unique(heatmap_data$ctgov)

# Heatmap of best performing models
heatmap_models <- heatmap_data %>% 
  ggplot(
    aes(
      x = ordered(ctgov, var_order_trials),
      y = fct_rev(cause),
      fill = dist
    )
  ) +
  geom_tile(colour = "black") +
  scale_fill_brewer(
    palette = "Accent",
    labels = c(
      "Exponential", "Gompertz", "Log-logistic",
      "Log-normal", "Weibull", "No events for cause"
    ),
    na.value = "white"
  ) +
  labs(
    x = NULL,
    y = NULL,
    fill = "Fitted distributions"
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(
      size = 12,
      angle = 45,
      hjust = 1,
      colour = "white"
    ),
    axis.text.y = element_text(
      size = 12,
      colour = "black"
    ),
    legend.title = element_text(
      size = 12,
      face = "bold",
      colour = "black"
    ),
    legend.text = element_text(
      size = 12,
      colour = "black"
    )
  )

# Save ggplot (with x-axis)
ggsave(
  "../Analysis/Output/plots/heatmap_models_coa_xaxis.pdf",
  heatmap_models,
  width = 7000,
  height = 1350,
  dpi = 300,
  units = "px"
)

# Save ggplot (w/o x-axis)
ggsave(
  "../Analysis/Output/plots/heatmap_models_coa.pdf",
  heatmap_models,
  width = 3500,
  height = 1350,
  dpi = 300,
  units = "px"
)

## Comparing best and second-best model fit per cause of attrition -------------

# Get table comparing AIC betwen best and second-best models
compare_aic <- import_fit %>%
  inner_join(import_summary %>% dplyr::select(ctgov)) %>% 
  dplyr::select(ctgov, cause, dist, aic) |> 
  filter(dist != "gengamma") |> 
  group_by(ctgov, cause) |> 
  slice_min(aic, n = 2) |> 
  slice_head(n = 2) |> 
  mutate(aic_diff = diff(aic)) |> 
  ungroup()

# Summarise AIC difference
compare_aic |> 
  dplyr::select(ctgov, cause, aic_diff) |> 
  distinct() |> 
  group_by(cause) |> 
  reframe(
    median = median(aic_diff),
    min = min(aic_diff),
    max = max(aic_diff)
  )

# Get a random 20% sample of trials per cause of attrition
# sample_trials_coa <- compare_aic |> 
#   dplyr::select(ctgov, cause) |> 
#   distinct() |> 
#   group_by(cause) |> 
#   sample_frac(0.2) |> 
#   ungroup()

# Join parameters to selected best and second-best models
join_aic_parameters <- compare_aic |> 
  inner_join(
    compare_aic |> 
      dplyr::select(ctgov, cause, dist)
  ) |> 
  inner_join(import_parameters) |> 
  dplyr::select(-se)

# Get KM estimates for these trials
km_trials <- compare_aic |> 
  dplyr::select(ctgov, cause) |> 
  distinct() |> 
  inner_join(import_t2e) |> 
  group_by(ctgov, cause) %>% 
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
  dplyr::select(ctgov, cause, time) |> 
  group_by(ctgov, cause) |> 
  summarise(
    min_time = min(time),
    max_time = max(time)
  ) |> 
  mutate(
    time = map2(min_time, max_time, seq)
  ) |> 
  unnest(cols = time) |> 
  ungroup() |> 
  dplyr::select(ctgov, cause, time) |> 
  left_join(km_trials)

# Estimate cdf at each time point (Exponential)
estimate_cdf_exp <- sequence_km_times |> 
  inner_join(
    join_aic_parameters |> 
      filter(dist == "exp") |> 
      dplyr::select(ctgov, cause, parameter, est) |> 
      pivot_wider(
        names_from = "parameter",
        values_from = "est"
      )
  ) |> 
  group_by(ctgov, time) |> 
  mutate(
    mdl_cdf = pexp(q = time, rate = exp(rate))
  ) |> 
  ungroup() |> 
  dplyr::select(-rate) |> 
  rename(
    km = estimate,
    exp = mdl_cdf
  ) |> 
  pivot_longer(
    km:exp,
    names_to = "type",
    values_to = "estimate"
  ) |> 
  na.omit()

# Estimate cdf at each time point (Gompertz)
estimate_cdf_gompertz <- sequence_km_times |> 
  inner_join(
    join_aic_parameters |> 
      filter(dist == "gompertz") |> 
      dplyr::select(ctgov, cause, parameter, est) |> 
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

# Estimate cdf at each time point (log-normal)
estimate_cdf_lnorm <- sequence_km_times |> 
  inner_join(
    join_aic_parameters |> 
      filter(dist == "lnorm") |> 
      dplyr::select(ctgov, cause, parameter, est) |> 
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

# Estimate cdf at each time point (log-logistic)
estimate_cdf_llogis <- sequence_km_times |> 
  inner_join(
    join_aic_parameters |> 
      filter(dist == "llogis") |> 
      dplyr::select(ctgov, cause, parameter, est) |> 
      pivot_wider(
        names_from = "parameter",
        values_from = "est"
      )
  ) |> 
  group_by(ctgov, time) |> 
  mutate(
    mdl_cdf = pllogis(q = time, scale = exp(scale), shape = exp(shape))
  ) |> 
  ungroup() |> 
  dplyr::select(-scale, -shape) |> 
  rename(
    km = estimate,
    llogis = mdl_cdf
  ) |> 
  pivot_longer(
    km:llogis,
    names_to = "type",
    values_to = "estimate"
  ) |> 
  na.omit()

# Estimate cdf at each time point (Weibull)
estimate_cdf_weibull <- sequence_km_times |> 
  inner_join(
    join_aic_parameters |> 
      filter(dist == "weibull") |> 
      dplyr::select(ctgov, cause, parameter, est) |> 
      pivot_wider(
        names_from = "parameter",
        values_from = "est"
      )
  ) |> 
  group_by(ctgov, time) |> 
  mutate(
    mdl_cdf = pweibull(q = time, shape = exp(shape), scale = exp(scale))
  ) |> 
  ungroup() |> 
  dplyr::select(-shape, -scale) |> 
  rename(
    km = estimate,
    weibull = mdl_cdf
  ) |> 
  pivot_longer(
    km:weibull,
    names_to = "type",
    values_to = "estimate"
  ) |> 
  na.omit()

# Combine both into one dataframe
estimate_cdf <- estimate_cdf_exp |> 
  full_join(estimate_cdf_gompertz) |> 
  full_join(estimate_cdf_llogis) |> 
  full_join(estimate_cdf_lnorm) |> 
  full_join(estimate_cdf_weibull) |> 
  arrange(ctgov, cause, time) |> 
  pivot_wider(
    names_from = "type",
    values_from = "estimate"
  ) %>% 
  pivot_longer(
    llogis:weibull,
    names_to = "dist",
    values_to = "model_est"
  )
  # mutate(
  #   new_id = dense_rank(ctgov),
  #   new_id = if_else(new_id < 10L, paste0("0", new_id), as.character(new_id)),
  #   new_id = paste0("Trial ", new_id),
  #   new_id = factor(new_id)
  # )

# Split into list by cause
lst_estimate_cdf <- split(estimate_cdf, estimate_cdf$cause)

# Vector for assigning colours in scale_colour_manual
vector_col <- c(
  "exp" = "skyblue",
  "gompertz" = "orange",
  "lnorm" = "purple",
  "llogis" = "deeppink",
  "weibull" = "red"
)

# Loop over list and plot
lst_plots <- list()
for (i in seq_along(lst_estimate_cdf)) {
  name <- names(lst_estimate_cdf[i])
  lst_plots[[name]] <- lst_estimate_cdf[[i]] %>% 
    filter(!is.na(model_est) & !is.na(km)) %>% 
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
        y = model_est,
        colour = dist
      ),
      linewidth = 0.5
    ) +
    scale_colour_manual(values = vector_col) +
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
      legend.position = "none"
    )
}

names(lst_plots) <- c("ae", "loe", "l2fu", "other", "pi", "pv", "vw")

# Save plots
for (i in seq_along(lst_plots)) {
  filepath <- "../Analysis/Output/plots/coa_compare_best_second_"
  name <- names(lst_plots[i])
  full_filepath <- paste0(filepath, name, ".png")
  print(name)
  ggsave(
    full_filepath,
    lst_plots[[i]],
    dpi = 300,
    width = 10,
    height = 6
  )
}

## Introductory statistics -----------------------------------------------------

# Conditions
import_flow %>% 
  group_by(condition) %>% 
  reframe(
    n = n(), 
    pcnt = n/90
  ) %>% 
  arrange(desc(pcnt)) %>% 
  write_csv("../Analysis/Output/summaries/aggregate_conditions.csv")

# Overall attrition statistics
import_summary %>% 
  dplyr::select(participants, dropouts) %>% 
  mutate(
    attrition = (dropouts/participants)*100
  ) %>% 
  reframe(
    total_attr_n = sum(dropouts),
    total_attr = round((sum(dropouts)/sum(participants)) * 100, 1),
    median_attr = round(median(attrition, na.rm = TRUE), 1),
    min_attr = round(min(attrition, na.rm = TRUE), 1),
    max_attr = round(max(attrition, na.rm = TRUE), 1),
    median_participants = median(participants, na.rm = TRUE),
    min_participants = min(participants, na.rm = TRUE),
    max_participants = max(participants, na.rm = TRUE)
  ) %>% 
  pivot_longer(
    everything(),
    names_to = "stat",
    values_to = "value"
  ) %>% 
  write_csv("../Analysis/Output/summaries/aggregate_flow.csv")

# Overall cause of attrition statistics
import_summary %>% 
  dplyr::select(dropouts:`Voluntary Withdrawal`) %>% 
  summarise(
    dropouts = sum(dropouts),
    across(
      `Adverse Event`:`Voluntary Withdrawal`,
      ~ sum(.)
    )
  ) %>% 
  pivot_longer(
    `Adverse Event`:`Voluntary Withdrawal`,
    names_to = "cause",
    values_to = "n"
  ) %>% 
  mutate(
    prop = round((n/dropouts) * 100, 1),
    value = paste0(n, " (", prop, ")")
  ) %>%
  dplyr::select(-dropouts) %>% 
  write_csv("../Analysis/Output/Summaries/aggregate_causes.csv")
  

# Table 1: Characteristics of included trials
import_summary %>% 
  dplyr::select(ctgov, condition, participants, dropouts) %>% 
  group_by(ctgov, condition) %>% 
  mutate(
    attrition = (dropouts/participants)*100
  ) %>% 
  group_by(condition) %>% 
  reframe(
    median_attr = round(median(attrition, na.rm = TRUE), 1),
    min_attr = round(min(attrition, na.rm = TRUE), 1),
    max_attr = round(max(attrition, na.rm = TRUE), 1),
    median_participants = median(participants, na.rm = TRUE),
    min_participants = min(participants, na.rm = TRUE),
    max_participants = max(participants, na.rm = TRUE)
  ) %>% 
  mutate(
    attr_value = paste0(median_attr, " (", min_attr, ", ", max_attr, ")"),
    prtcp_value = paste0(median_participants, " (", min_participants, ", ", max_participants, ")")
  ) %>% 
  write_csv("../Analysis/Output/Summaries/summary_tbl1_total.csv")

# Table 2: Frequency of cause-specific attrition among included trials
import_summary %>% 
  dplyr::select(condition, dropouts:`Voluntary Withdrawal`) %>% 
  pivot_longer(
    `Adverse Event`:`Voluntary Withdrawal`,
    names_to = "cause",
    values_to = "n"
  ) %>% 
  group_by(condition) %>% 
  mutate(
    prop = (n/dropouts)*100
  ) %>% 
  group_by(condition, cause) %>% 
  reframe(
    median = round(median(prop, na.rm = TRUE), 1),
    min = round(min(prop, na.rm = TRUE), 1),
    max = round(max(prop, na.rm = TRUE), 1)
  ) %>% 
  mutate(
    value = paste0(median, " (", min, ", ", max, ")")
  ) %>% 
  dplyr::select(-median, -min, -max) %>% 
  pivot_wider(
    names_from = "cause",
    values_from = "value"
  ) %>% 
  write_csv("../Analysis/Output/Summaries/summary_tbl1_causes.csv")

# Supplementary Table 2: Frequency of causes of attrition among individual trials
import_summary %>% 
  dplyr::select(ctgov, condition, dropouts:`Voluntary Withdrawal`) %>% 
  pivot_longer(
    `Adverse Event`:`Voluntary Withdrawal`,
    names_to = "cause",
    values_to = "n"
  ) %>% 
  group_by(ctgov) %>% 
  mutate(
    pcnt = round((n/dropouts) * 100, 1),
    count = paste0(n, "/", dropouts, " (", pcnt, ")")
  ) %>% 
  dplyr::select(ctgov, condition, cause, count) %>% 
  pivot_wider(
    names_from = "cause",
    values_from = "count"
  ) %>% 
  ungroup() %>% 
  write_csv("../Analysis/Output/Summaries/summary_suppltbl2.csv")

# Model outputs for github
github_part1 <- unnest(data$best_models, c("nested_params")) %>% 
  filter(
    !(ctgov %in% c("NCT01131676", "NCT00274573"))
  ) %>% 
  inner_join(
    data$parameters_fit %>% 
      dplyr::select(ctgov, cause, dist, aic)
  )

# Unlist parameters and join to model outputs
lst_params <- github_part1$params

for (i in seq_along(github_part1$params)){
  lst_params[[i]] <- as.data.frame(lst_params[[i]])
}

df_params <- bind_rows(lst_params)

# Pivot to long data
model_parameters <- github_part1 %>% 
  cbind(df_params) %>% 
  pivot_longer(
    meanlog:scale,
    names_to = "parameter",
    values_to = "mean"
  ) %>% 
  filter(!is.na(mean)) %>% 
  dplyr::select(-params) %>% 
  inner_join(
    import_summary %>% 
      dplyr::select(ctgov, condition)
  ) %>% 
  arrange(condition, ctgov, cause) %>% 
  dplyr::select(ctgov, condition, cause, dist, aic, parameter, mean) %>% 
  mutate(
    dist = case_when(
      dist == "exp" ~ "Exponential",
      dist == "gengamma" ~ "Generalised gamma",
      dist == "gompertz" ~ "Gompertz",
      dist == "lnorm" ~ "Log-normal",
      dist == "weibull" ~ "Weibull"
    ),
    parameter = case_when(
      parameter == "rate" ~ "Rate",
      parameter == "shape" ~ "Shape",
      parameter == "meanlog" ~ "Mean log",
      parameter == "sdlog" ~ "SD log",
      parameter == "scale" ~ "Scale",
      TRUE ~ parameter
    ),
    mean = round(mean, 3)
  )

write_csv(model_parameters, "../Analysis/Output/github/model_outputs.csv")

## Summarise hazard rates ------------------------------------------------------

# Import tally
import_tally_hazards <- read_csv("../Analysis/Output/coa_tbls/tally_hazards_causes.csv") %>% 
  rename(var = `ctgov\tcause\tpeak\tdirection\tshape`) %>% 
  separate(
    var,
    into = c("ctgov", "cause", "peak", "direction", "shape"), sep = "\t"
  ) %>% 
  inner_join(import_summary %>% dplyr::select(ctgov))

# High point
tbl_hazard_highest <- import_tally_hazards %>% 
  group_by(cause, dist, peak) %>% 
  reframe(
    n = n()
  ) %>% 
  group_by(cause) %>% 
  mutate(
    cause_n = case_when(
      cause == "Adverse Event" ~ 81,
      cause == "Lack of Efficacy" ~ 20,
      cause == "Lost to Follow-up" ~ 24,
      cause == "Other/Miscellaneous" ~ 40,
      cause == "PI/Sponsor Decision" ~ 6,
      cause == "Protocol Violation" ~ 28,
      cause == "Voluntary Withdrawal" ~ 64
    ),
    pcnt = round((n/cause_n) * 100, 1),
    value = paste0(n, " (", pcnt, ")")
  ) %>% 
  group_by(cause, peak) %>% 
  mutate(
    total = sum(n),
    total_pcnt = round(total/cause_n * 100, 1),
    total_value = paste0(total, " (", total_pcnt, ")")
  ) %>% 
  ungroup() %>% 
  pivot_wider(
    names_from = "peak",
    values_from = "value"
  ) %>% 
  write_csv(
    "../Analysis/Output/coa_tbls/hazards_peak.csv"
  )

# Model and direction
tbl_hazard_model_direction <- import_tally_hazards %>% 
  group_by(cause, dist, direction) %>% 
  reframe(
    n = n()
  ) %>% 
  group_by(cause) %>% 
  mutate(
    cause_n = case_when(
      cause == "Adverse Event" ~ 81,
      cause == "Lack of Efficacy" ~ 20,
      cause == "Lost to Follow-up" ~ 24,
      cause == "Other/Miscellaneous" ~ 40,
      cause == "PI/Sponsor Decision" ~ 6,
      cause == "Protocol Violation" ~ 28,
      cause == "Voluntary Withdrawal" ~ 64
    ),
    pcnt = round((n/cause_n) * 100, 1),
    value = paste0(n, " (", pcnt, ")")
  ) %>%
  group_by(cause, direction) %>% 
  mutate(
    total = sum(n),
    total_pcnt = round(total/cause_n * 100, 1),
    total_value = paste0(total, " (", total_pcnt, ")")
  ) %>% 
  ungroup() %>% 
  pivot_wider(
    names_from = "direction",
    values_from = "value"
  ) %>% 
  write_csv(
    "../Analysis/Output/coa_tbls/hazards_direction.csv"
  )

# Model and shape
tbl_hazard_model_shape <- import_tally_hazards %>% 
  group_by(cause, dist, shape) %>% 
  reframe(
    n = n()
  ) %>% 
  group_by(cause) %>% 
  mutate(
    cause_n = case_when(
      cause == "Adverse Event" ~ 81,
      cause == "Lack of Efficacy" ~ 20,
      cause == "Lost to Follow-up" ~ 24,
      cause == "Other/Miscellaneous" ~ 40,
      cause == "PI/Sponsor Decision" ~ 6,
      cause == "Protocol Violation" ~ 28,
      cause == "Voluntary Withdrawal" ~ 64
    ),
    pcnt = round((n/cause_n) * 100, 1),
    value = paste0(n, " (", pcnt, ")")
  ) %>%
  group_by(cause, shape) %>% 
  mutate(
    total = sum(n),
    total_pcnt = round(total/cause_n * 100, 1),
    total_value = paste0(total, " (", total_pcnt, ")")
  ) %>% 
  ungroup() %>% 
  pivot_wider(
    names_from = "shape",
    values_from = "value"
  ) %>% 
  write_csv(
    "../Analysis/Output/coa_tbls/hazards_shape.csv"
  )


