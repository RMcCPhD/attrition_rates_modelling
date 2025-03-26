
# Plot estimates from empirical data and models
# Plot files are not saved in repo to save file space, but can be obtained
# by running through the scripts

source("config.R")
library(tidyverse)
library(ggh4x)
library(forcats)
library(flexsurv)

# Import data
imp_params <- read_csv("total_attrition/data/parameters.csv")
lst_prep <- read_rds("total_attrition/processed_data/prepared_data.rds")
imp_cind <- read_rds("total_attrition/processed_data/cind_sampled_df.rds")
imp_haz <- read_rds("total_attrition/processed_data/haz_sampled_df.rds")
imp_sum <- read_csv("vivli/data/summaries/participant_flow_formatted.csv")

# Create abbreviations for conditions
cond_short <- imp_sum %>% 
  group_by(condition) %>% 
  reframe() %>% 
  mutate(
    condition_short = case_match(
      condition,
      "Asthma" ~ "Asth",
      "Benign prostatic hyperplasia" ~ "BPH",
      "Hypertension" ~ "Hyp",
      "Osteoarthritis" ~ "Oarth",
      "Osteoporosis" ~ "Oporo",
      "Parkinson's disease" ~ "Park",
      "Pulmonary fibrosis" ~ "Pfib",
      "Restless leg syndrome" ~ "RLS",
      "Type 2 diabetes" ~ "T2DM",
      .default = condition
    )
  )

# Create index for grouping by condition and sample size
ind_cond_n <- imp_sum %>% 
  select(ctgov, participants, condition) %>% 
  rename(n = participants) %>% 
  left_join(cond_short) %>% 
  arrange(condition, desc(n)) %>% 
  group_by(condition) %>% 
  mutate(
    trial_num = row_number(),
    trial_num = if_else(trial_num < 10, paste0("0", trial_num), as.character(trial_num)),
    plot_id = paste0(condition_short, " ", trial_num, ", n=", n)
  ) %>% 
  ungroup() %>% 
  select(-trial_num)

# Join new id to cumulative incidence and hazard rate datasets
# Remove a few iterations that make a plot unreadable (COPD_19_n=472), has
# 96 sampled trajectories
cind_plot_df <- ind_cond_n %>% 
  left_join(imp_cind) %>% 
  filter(!is.na(estimate)) %>% 
  filter(!(plot_id == "COPD 19, n=472" & iter %in% c(c(16L, 39L, 49L, 56L))))

haz_plot_df <- ind_cond_n %>% 
  left_join(imp_haz) %>% 
  filter(dist != "exp")

### Heatmap ####################################################################

# Figure 1 (heatmap of best model performance compared to others)

# Get data for heatmap
imp_hmap_data <- read_csv("total_attrition/data/parameters.csv") %>% 
  select(ctgov, dist) %>% 
  distinct() %>% 
  filter(!(ctgov %in% c("NCT01131676", "NCT00274573"))) %>% 
  left_join(read_csv("total_attrition/data/fit.csv")) %>% 
  select(-loglik) %>% 
  arrange(ctgov)

# Identify best performing model as centre for aic scale
imp_hmap_id_best <- imp_hmap_data %>% 
  left_join(
    read_csv("usable_data/total_attrition/params_pblc.csv") %>% 
      select(ctgov, best = dist) %>% 
      distinct()
  ) %>% 
  filter(dist != "gengamma") %>% 
  mutate(flag = if_else(best == dist, 1, 0)) %>% 
  select(-best) %>% 
  group_by(ctgov) %>% 
  mutate(
    best_aic = aic[flag == 1],
    max_aic = max(aic),
    min_aic = min(aic),
    aic_scale = pmin(aic - min_aic, 10)
  ) %>% 
  ungroup() %>% 
  mutate(
    aic_scale_colour = case_when(
      aic_scale == 0 ~ "dark grey", # Best model(s)
      aic_scale <= 2 ~ "green", # AIC diff of <= 2 from best
      aic_scale <= 4 ~ "yellow", # AIC diff of <= 4 from best
      aic_scale <= 7 ~ "orange", # AIC diff of <= 7 from best
      TRUE ~ "white" # Worst (AIC diff >7 from best)
    )
  )

# Prepare heatmap data
hmap_df <- ind_cond_n %>% 
  left_join(imp_hmap_id_best) %>% 
  select(plot_id, condition, dist, aic_scale, aic_scale_colour) %>% 
  group_by(condition) %>% 
  mutate(
    plot_index = as.integer(factor(plot_id)),
  ) %>% 
  ungroup() %>% 
  mutate(
    dist = case_match(
      dist,
      "weibull" ~ "Weibull",
      "lnorm" ~ "Log-normal",
      "llogis" ~ "Log-logistic",
      "gompertz" ~ "Gompertz",
      "exp" ~ "Exponential"
    ),
    dist = factor(
      dist,
      levels = c("Gompertz", "Log-normal", "Log-logistic", "Weibull", "Exponential")
    )
  )

# Heatmap
hmap <- hmap_df %>% 
  ggplot(aes(x = dist, y = fct_rev(plot_id), fill = aic_scale_colour)) +
  geom_tile(colour = "black") +
  scale_fill_identity() +
  labs(x = NULL, y = NULL) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 8, angle = 45, hjust = 1, colour = "black"),
    axis.text.y = element_text(colour = "white", size = 10)
  )

# ggsave(
#   "total_attrition/plots/fig1_raw.png",
#   hmap,
#   dpi = 300,
#   width = 90,
#   height = 120,
#   units = "mm"
# )

# ggsave(
#   "total_attrition/plots/fig1_raw_labels.png",
#   hmap,
#   dpi = 300,
#   width = 90,
#   height = 360,
#   units = "mm"
# )

### Cumulative incidence #######################################################

# Figure 2 (empirical vs fitted estimates)
plot_cind <- cind_plot_df %>% 
  ggplot(aes(x = time, y = round((1 - estimate) * 100, 1), group = iter)) +
  geom_step(linewidth = 0.5, colour = "black") +
  geom_line(
    aes(time, round((1 - fit_est) * 100, 1)), 
    linewidth = 0.5, 
    alpha = 0.02,
    colour = case_when(
      cind_plot_df$dist == "exp" ~ "blue",
      cind_plot_df$dist == "gompertz" ~"orange",
      cind_plot_df$dist == "lnorm" ~ "purple",
      cind_plot_df$dist == "llogis" ~ "deeppink",
      cind_plot_df$dist == "weibull" ~ "red"
    )
  ) +
  labs(x = "Time in days", y = "Cumulative incidence of attrition") +
  facet_wrap(~plot_id, scales = "free") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
  theme_bw() +
  theme(text = element_text(size = 5))

# ggsave(
#   "total_attrition/plots/fig2.png",
#   plot_cind,
#   dpi = 300,
#   width = 190,
#   height = 160,
#   units = "mm"
# )

## Supplementary figure 3 (confidence limits)
# Get output parameters and standard errors
ci_params <- lst_prep$params %>% 
  unnest(cols = "params") %>% 
  mutate(params = map(params, as_tibble)) %>% 
  unnest(params) %>% 
  select(ctgov, dist) %>% 
  left_join(imp_params) %>% 
  pivot_wider(names_from = "parameter", values_from = c("est", "se"))

# Sequence all time points per trial and join empirical estimates
ci_cind <- lst_prep$cind %>% 
  unnest(cols = "cind") %>% 
  group_by(ctgov) %>% 
  mutate(time_seq = list(seq(min(time), max(time), by = 1))) %>% 
  unnest(time_seq) %>% 
  select(ctgov, time_seq) %>% 
  left_join(
    lst_prep$cind %>% 
      unnest(cols = "cind") %>% 
      rename(time_seq = time)
  ) %>% 
  distinct()

# Estimate model cumulative incidence and confidence limits
ci_cind_mdls <- ci_cind %>% 
  left_join(ci_params, by = "ctgov") %>% 
  group_by(ctgov) %>% 
  mutate(
    mdl = pmap_dbl(
      list(time_seq, dist, est_rate, est_shape, est_scale, est_meanlog, est_sdlog),
      ~ {
        time <- ..1
        dist <- ..2
        rate <- ..3
        shape <- ..4
        scale <- ..5
        meanlog <- ..6
        sdlog <- ..7
        
        if (dist == "exp") {
          1 - pexp(time, rate = exp(rate))
          
        } else if (dist == "gompertz") {
          1 - pgompertz(time, shape = shape, rate = exp(rate))
          
        } else if (dist == "llogis") {
          1 - pllogis(time, shape = exp(shape), scale = exp(scale))
          
        } else if (dist == "lnorm") {
          1 - plnorm(time, meanlog = meanlog, sdlog = exp(sdlog))
          
        } else if (dist == "weibull") {
          1 - pweibull(time, shape = exp(shape), scale = exp(scale))
          
        } else {
          NA_real_
        }
      }
    ),
    se_mdl = pmap_dbl(
      list(
        time_seq, dist, est_rate, est_shape, est_scale, est_meanlog, est_sdlog, 
        se_rate, se_shape, se_scale, se_meanlog, se_sdlog
      ),
      ~ {
        time <- ..1
        dist <- ..2
        rate <- ..3
        shape <- ..4
        scale <- ..5
        meanlog <- ..6
        sdlog <- ..7
        se_rate <- ..8
        se_shape <- ..9
        se_scale <- ..10
        se_meanlog <- ..11
        se_sdlog <- ..12
        
        if (dist == "exp") {
          S_t <- exp(-exp(rate) * time)
          grad_rate <- -time * S_t * exp(rate)
          se_S_t <- sqrt(grad_rate^2 * se_rate^2)
          
        } else if (dist == "gompertz") {
          S_t <- exp(-(exp(rate) / shape) * (exp(shape * time) - 1))
          grad_rate <- -S_t * (1 / shape) * (exp(shape * time) - 1) * exp(rate)
          grad_shape <- S_t * (exp(rate) / shape^2) * (exp(shape * time) - 1) - S_t * (exp(rate) / shape) * exp(shape * time) * time
          se_S_t <- sqrt(grad_rate^2 * se_rate^2 + grad_shape^2 * se_shape^2)
          
        } else if (dist == "llogis") {
          S_t <- 1 / (1 + (time / exp(scale))^exp(shape))
          grad_shape <- -S_t^2 * log(time / exp(scale)) * (time / exp(scale))^exp(shape) * exp(shape)
          grad_scale <- S_t^2 * exp(shape) * (time / exp(scale))^exp(shape) / exp(scale)
          se_S_t <- sqrt(grad_shape^2 * se_shape^2 + grad_scale^2 * se_scale^2)
          
        } else if (dist == "lnorm") {
          S_t <- 1 - plnorm(time, meanlog, exp(sdlog))
          grad_meanlog <- -dnorm((log(time) - meanlog) / exp(sdlog)) / exp(sdlog)
          grad_sdlog <- -dnorm((log(time) - meanlog) / exp(sdlog)) * (log(time) - meanlog) / (exp(sdlog)^2) * exp(sdlog)
          se_S_t <- sqrt(grad_meanlog^2 * se_meanlog^2 + grad_sdlog^2 * se_sdlog^2)
          
        } else if (dist == "weibull") {
          S_t <- exp(-(time / exp(scale))^exp(shape))
          grad_shape <- -S_t * log(time / exp(scale)) * (time / exp(scale))^exp(shape) * exp(shape)
          grad_scale <- S_t * exp(shape) * (time / exp(scale))^exp(shape) / exp(scale)
          se_S_t <- sqrt(grad_shape^2 * se_shape^2 + grad_scale^2 * se_scale^2)
          
        } else {
          NA_real_
        }
      }
    ),
    
    # Compute CI based on delta method
    lci = pmax((mdl + 1.96 * se_mdl), 0),
    uci = pmin((mdl - 1.96 * se_mdl), 1)
  ) %>% 
  select(ctgov:dist, mdl:uci)

# Add plot ids
ci_cind_plot_df <- ind_cond_n %>% 
  left_join(ci_cind_mdls) %>% 
  filter(!is.na(estimate))

# Plot
plot_cind_ci <- ci_cind_plot_df %>% 
  ggplot(aes(x = time_seq, y = round((1 - estimate) * 100, 1))) +
  geom_step(linewidth = 0.5, colour = "black") +
  geom_line(
    aes(x = time_seq, y = round((1 - mdl) * 100, 1)),
    linewidth = 0.5,
    colour = case_when(
      ci_cind_plot_df$dist == "exp" ~ "blue",
      ci_cind_plot_df$dist == "gompertz" ~"orange",
      ci_cind_plot_df$dist == "lnorm" ~ "purple",
      ci_cind_plot_df$dist == "llogis" ~ "deeppink",
      ci_cind_plot_df$dist == "weibull" ~ "red"
    )
  ) +
  geom_ribbon(
    aes(
      x = time_seq,
      ymin = round((1 - lci) * 100, 1),
      ymax = round((1 - uci) * 100, 1)
    ),
    alpha = 0.4,
    fill = case_when(
      ci_cind_plot_df$dist == "exp" ~ "skyblue",
      ci_cind_plot_df$dist == "gompertz" ~"orange",
      ci_cind_plot_df$dist == "lnorm" ~ "purple",
      ci_cind_plot_df$dist == "llogis" ~ "deeppink",
      ci_cind_plot_df$dist == "weibull" ~ "red"
    )
  ) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
  labs(x = "Time in days", y = "Cumulative incidence of attrition") +
  facet_wrap(~plot_id, scales = "free") +
  theme_bw() +
  theme(text = element_text(size = 6))

ggsave(
  "total_attrition/plots/suppfig2.png",
  plot_cind_ci,
  dpi = 300,
  width = 10,
  height = 6
)

## Supplementary figure 3 (Gompertz vs log-normal)
# Get Gompertz and log-normal parameters where either were best-fitting
# 57 trials
gmp_lnorm_params <- ci_cind_mdls %>% 
  count(ctgov, dist) %>% 
  filter(dist %in% c("gompertz", "lnorm")) %>% 
  select(ctgov) %>% 
  left_join(imp_params %>% filter(dist %in% c("gompertz", "lnorm"))) %>% 
  select(-se) %>% 
  pivot_wider(names_from = c("dist", "parameter"), values_from = c("est"))

gmp_lnorm_mdls <- ci_cind %>% 
  filter(ctgov %in% gmp_lnorm_params$ctgov) %>% 
  left_join(gmp_lnorm_params, by = "ctgov") %>% 
  group_by(ctgov) %>% 
  mutate(
    mdl_gmp = 1 - pgompertz(time_seq, shape = gompertz_shape, rate = exp(gompertz_rate)),
    mdl_lnorm = 1 - plnorm(time_seq, meanlog = lnorm_meanlog, sdlog = exp(lnorm_sdlog))
  )

plot_gmp_lnorm <- ind_cond_n %>% 
  filter(ctgov %in% gmp_lnorm_mdls$ctgov) %>% 
  left_join(gmp_lnorm_mdls) %>% 
  filter(!is.na(estimate)) %>% 
  ggplot(aes(x = time_seq, y = 1 - estimate)) +
  geom_step(colour = "black", linewidth = 0.5) +
  geom_line(aes(x = time_seq, y = 1 - mdl_gmp), colour = "orange", linewidth = 0.5) +
  geom_line(aes(x = time_seq, y = 1 - mdl_lnorm), colour = "purple", linewidth = 0.5) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
  labs(x = "Time in days", y = "Cumulative incidence of attrition") +
  facet_wrap(~plot_id, scales = "free", ncol = 9) +
  theme_bw() +
  theme(text = element_text(size = 6))

ggsave(
  "total_attrition/plots/suppfig3.png",
  plot_gmp_lnorm,
  dpi = 300,
  width = 10,
  height = 6
)


### Hazard rate ################################################################

# Figure 3
plot_haz <- haz_plot_df %>% 
  group_by(ctgov) %>% 
  mutate(est = haz_est/max(haz_est)) %>% 
  ggplot(aes(x = time, y = est)) +
  geom_line(
    linewidth = 0.5,
    colour = case_when(
      haz_plot_df$dist == "gompertz" ~"orange",
      haz_plot_df$dist == "lnorm" ~ "purple",
      haz_plot_df$dist == "llogis" ~ "deeppink",
      haz_plot_df$dist == "weibull" ~ "red"
    )
  ) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
  labs(x = "Time in days", y = "Scaled hazard rate") +
  facet_wrap(~plot_id, scales = "free") +
  theme_bw() +
  theme(text = element_text(size = 5))

ggsave(
  "total_attrition/plots/fig3.png",
  plot_haz,
  dpi = 300,
  width = 190,
  height = 160,
  units = "mm"
)
