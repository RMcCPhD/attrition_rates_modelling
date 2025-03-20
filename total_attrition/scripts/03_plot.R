
# Plot estimates from empirical data and models
# Plot files are not saved in repo to save file space, but can be obtained
# by running through the scripts

source("config.R")
library(tidyverse)
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
    plot_id = paste0(condition_short, "_", trial_num, "_n=", n)
  ) %>% 
  select(-trial_num)

# Address error for NCT01694771, randomised period ends at 84 days but follow-up
# period was a threshold for completion (i.e. ~22 days), actual end point is
# 84 days
imp_cind_fix <- as_tibble(imp_cind) %>% 
  group_by(ctgov) %>% 
  mutate(
    estimate = case_when(
      ctgov == "NCT01694771" 
      & time > 83
      & !is.na(estimate) ~ first(estimate[time == 83]),
      TRUE ~ estimate
    )
  ) %>% 
  ungroup()

# Join new id to cumulative incidence and hazard rate datasets
cind_plot_df <- ind_cond_n %>% 
  left_join(imp_cind_fix) %>% 
  filter(!is.na(estimate))

haz_plot_df <- ind_cond_n %>% 
  left_join(imp_haz) %>% 
  filter(dist != "exp")

### Cumulative incidence #######################################################

# Figure 1 (empirical vs fitted estimates)
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
  scale_x_continuous(n.breaks = 4) +
  scale_y_continuous(n.breaks = 4) +
  labs(x = "Time in days", y = "Cumulative incidence of attrition") +
  facet_wrap(~plot_id, scales = "free") +
  theme_bw() +
  theme(text = element_text(size = 5))

ggsave(
  "total_attrition/plots/fig1.png",
  plot_cind,
  dpi = 300,
  width = 2244/300,
  height = 1683/300,
  units = "in"
)

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
      group_by(ctgov) %>% 
      mutate(
        estimate = case_when(
          ctgov == "NCT01694771" 
          & time > 83
          & !is.na(estimate) ~ first(estimate[time == 83]),
          TRUE ~ estimate
        )
      ) %>% 
      ungroup() %>% 
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
ci_cind_plot_df <- ind_cond_n %>% left_join(ci_cind_mdls) %>% filter(!is.na(estimate))

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
  scale_x_continuous(n.breaks = 4) +
  scale_y_continuous(n.breaks = 4) +
  labs(x = "Time in days", y = "Cumulative incidence of attrition") +
  facet_wrap(~plot_id, scales = "free") +
  theme_bw() +
  theme(text = element_text(size = 5))

ggsave(
  "total_attrition/plots/suppfig3.png",
  plot_cind_ci,
  dpi = 300,
  width = 2244/300,
  height = 1683/300,
  units = "in"
)

## Supplementary figure 4 (Gompertz vs log-normal)


### Hazard rate ################################################################

# Figure 2
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
  scale_x_continuous(n.breaks = 4) +
  scale_y_continuous(n.breaks = 4) +
  labs(x = "Time in days", y = "Scaled hazard rate") +
  facet_wrap(~plot_id, scales = "free") +
  theme_bw() +
  theme(text = element_text(size = 5))

ggsave(
  "total_attrition/plots/fig2.png",
  plot_haz,
  dpi = 300,
  width = 2244/300,
  height = 1683/300,
  units = "in"
)
