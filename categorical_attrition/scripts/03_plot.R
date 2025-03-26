
# Plot estimates from empirical data and models
# Plot files are not saved in repo to save file space, but can be obtained
# by running through the scripts

source("config.R")
library(tidyverse)
library(flexsurv)

# Import data
imp_params <- read_csv("categorical_attrition/data/parameters.csv")
lst_prep <- read_rds("categorical_attrition/processed_data/prepared_data.rds")
imp_cind <- read_rds("categorical_attrition/processed_data/cind_sampled_df.rds")
imp_haz <- read_rds("categorical_attrition/processed_data/haz_sampled_df.rds")
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
cind_plot_df <- ind_cond_n %>% left_join(imp_cind)
haz_plot_df <- ind_cond_n %>% left_join(imp_haz) %>% filter(dist != "exp")

### Heatmap ####################################################################

# Figure 1 (heatmap of best model performance compared to others)

# Get data for heatmap
imp_hmap_data <- read_csv("categorical_attrition/data/parameters.csv") %>% 
  select(ctgov, dist) %>% 
  distinct() %>% 
  filter(!(ctgov %in% c("NCT01131676", "NCT00274573"))) %>% 
  left_join(read_csv("categorical_attrition/data/fit.csv")) %>% 
  select(-loglik) %>% 
  arrange(ctgov)

# Identify best performing model as centre for aic scale
imp_hmap_id_best <- imp_hmap_data %>% 
  left_join(
    read_csv("usable_data/categorical_attrition/params_pblc.csv") %>% 
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
    aic_scale = (max_aic - aic) / (max_aic - min_aic)
  ) %>% 
  ungroup()

# Prepare heatmap data
hmap_df <- ind_cond_n %>% 
  left_join(imp_hmap_id_best) %>% 
  select(plot_id, condition, dist, aic_scale) %>% 
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
    )
  )

# Heatmap
hmap <- hmap_df %>% 
  ggplot(aes(x = plot_index, y = fct_rev(dist), fill = aic_scale)) +
  geom_tile(colour = "black") +
  scale_x_continuous(
    breaks = function(x) seq(floor(min(x)), ceiling(max(x)), by = 1),
    expand = expansion(mult = c(0, 0.01))
  ) +
  scale_fill_gradient(
    low = "white", 
    high = "steelblue",
    name = "Scaled AIC",
    breaks = c(0, 0.25, 0.75, 1),
    labels = c("0 (Worst)", "", "", "1 (Best)")
  ) +
  facet_wrap(~ condition, scales = "free_x", ncol = 3) +
  labs(x = NULL, y = NULL) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank()
  )

ggsave(
  "categorical_attrition/plots/fig1.png",
  hmap,
  dpi = 300,
  width = 2244/300,
  height = 1683/300,
  units = "in"
)

### Hazard rate ################################################################

# Figure 5
plot_haz <- haz_plot_df %>% 
  group_by(ctgov, cause) %>% 
  mutate(est = haz_est/max(haz_est)) %>% 
  ggplot(aes(x = time, y = est, group = plot_id, colour = dist)) +
  geom_line(linewidth = 0.5) +
  scale_colour_manual(
    values = c(
      "gompertz" = "orange",
      "lnorm" = "purple",
      "llogis" = "deeppink",
      "weibull" = "red"
    )
  ) +
  scale_x_continuous(n.breaks = 4) +
  scale_y_continuous(n.breaks = 4) +
  labs(x = "Time in days", y = "Scaled hazard rate", colour = NULL) +
  facet_wrap(~cause, scales = "free") +
  theme_bw() +
  theme(
    text = element_text(size = 6),
    legend.position = "none"
  )

ggsave(
  "categorical_attrition/plots/fig5.png",
  plot_haz,
  dpi = 300,
  width = 2244/300,
  height = 1683/300,
  units = "in"
)

### Cumulative incidence #######################################################

# Adverse events (81 trials)
plot_cind <- cind_plot_df %>% 
  filter(cause == "Adverse Event") %>% 
  mutate(across(c(estimate, fit_est), ~ round((1 - .) * 100, 1))) %>% 
  group_by(plot_id, iter) %>% 
  filter(
    plot_id != "BPH 03, n=606" | (plot_id == "BPH 03, n=606" & all(fit_est <= 3)),
    plot_id != "COPD 19, n=472" | (plot_id == "COPD 19, n=472" & all(fit_est <= 6)),
    plot_id != "Hyp 01, n=1039" | (plot_id == "Hyp 01, n=1039" & all(fit_est <= 4)),
    plot_id != "Hyp 08, n=426" | (plot_id == "Hyp 08, n=426" & all(fit_est <= 5)),
    plot_id != "Park 04, n=475" | (plot_id == "Park 04, n=475" & all(fit_est <= 7)),
    plot_id != "T2DM 04, n=1360" | (plot_id == "T2DM 04, n=1360" & all(fit_est <= 4)),
    plot_id != "T2DM 11, n=899" | (plot_id == "T2DM 11, n=899" & all(fit_est <= 5)),
    plot_id != "T2DM 22, n=574" | (plot_id == "T2DM 22, n=574" & all(fit_est <= 6)),
    plot_id != "T2DM 25, n=561" | (plot_id == "T2DM 25, n=561" & all(fit_est <= 7)),
    plot_id != "T2DM 28, n=492" | (plot_id == "T2DM 28, n=492" & all(fit_est <= 5)),
    plot_id != "T2DM 29, n=447" | (plot_id == "T2DM 29, n=447" & all(fit_est <= 7)),
    plot_id != "T2DM 30, n=389" | (plot_id == "T2DM 30, n=389" & all(fit_est <= 5)),
    plot_id != "T2DM 31, n=302" | (plot_id == "T2DM 31, n=302" & all(fit_est <= 7)),
    plot_id != "T2DM 34, n=299" | (plot_id == "T2DM 34, n=299" & all(fit_est <= 7))
  ) %>% 
  ungroup() %>% 
  ggplot(aes(x = time, group = iter)) +
  geom_step(
    data = ~ filter(.x, !is.na(estimate)), 
    aes(y = estimate), 
    linewidth = 0.5, 
    colour = "black"
  ) +
  geom_line(aes(y = fit_est, colour = dist), linewidth = 0.5, alpha = 0.02) +
  scale_colour_manual(
    values = c(
      "exp" = "blue",
      "gompertz" = "orange",
      "lnorm" = "purple",
      "llogis" = "deeppink",
      "weibull" = "red"
    )
  ) +
  scale_x_continuous(n.breaks = 4) +
  scale_y_continuous(n.breaks = 4) +
  labs(x = "Time in days", y = "Cumulative incidence of attrition") +
  facet_wrap(~plot_id, scales = "free") +
  theme_bw() +
  theme(
    text = element_text(size = 5),
    legend.position = "none"
  )

ggsave(
  "categorical_attrition/plots/suppl_ae.png",
  plot_cind,
  dpi = 300,
  width = 2244/300,
  height = 1683/300,
  units = "in"
)

# Lack of efficacy (20 trials)
plot_cind <- cind_plot_df %>% 
  filter(cause == "Lack of Efficacy") %>% 
  mutate(across(c(estimate, fit_est), ~ round((1 - .) * 100, 1))) %>% 
  group_by(plot_id, iter) %>%
  filter(
    plot_id != "BPH 03, n=606" | (plot_id == "BPH 03, n=606" & all(fit_est <= 6)),
    plot_id != "COPD 03, n=2488" | (plot_id == "COPD 03, n=2488" & all(fit_est <= 1)),
    plot_id != "COPD 14, n=624" | (plot_id == "COPD 14, n=624" & all(fit_est <= 6)),
    plot_id != "COPD 18, n=519" | (plot_id == "COPD 18, n=519" & all(fit_est <= 6)),
    plot_id != "Hyp 07, n=490" | (plot_id == "Hyp 07, n=490" & all(fit_est <= 10)),
    plot_id != "T2DM 23, n=566" | (plot_id == "T2DM 23, n=566" & all(fit_est <= 10))
  ) %>%
  ungroup() %>%
  ggplot(aes(x = time, group = iter)) +
  geom_step(
    data = ~ filter(.x, !is.na(estimate)), 
    aes(y = estimate), 
    linewidth = 0.5, 
    colour = "black"
  ) +
  geom_line(aes(y = fit_est, colour = dist), linewidth = 0.5, alpha = 0.02) +
  scale_colour_manual(
    values = c(
      "exp" = "blue",
      "gompertz" = "orange",
      "lnorm" = "purple",
      "llogis" = "deeppink",
      "weibull" = "red"
    )
  ) +
  scale_x_continuous(n.breaks = 4) +
  scale_y_continuous(n.breaks = 4) +
  labs(x = "Time in days", y = "Cumulative incidence of attrition") +
  facet_wrap(~plot_id, scales = "free") +
  theme_bw() +
  theme(
    text = element_text(size = 5),
    legend.position = "none"
  )

ggsave(
  "categorical_attrition/plots/suppl_loe.png",
  plot_cind,
  dpi = 300,
  width = 2244/300,
  height = 1683/300,
  units = "in"
)

# Lost to follow-up (24 trials)
plot_cind <- cind_plot_df %>% 
  filter(cause == "Lost to Follow-up") %>% 
  mutate(across(c(estimate, fit_est), ~ round((1 - .) * 100, 1))) %>% 
  group_by(plot_id, iter) %>%
  filter(
    plot_id != "COPD 04, n=1829" | (plot_id == "COPD 04, n=1829" & all(fit_est <= 2)),
    plot_id != "COPD 08, n=983" | (plot_id == "COPD 08, n=983" & all(fit_est <= 5)),
    plot_id != "Hyp 04, n=860" | (plot_id == "Hyp 04, n=860" & all(fit_est <= 6)),
    plot_id != "T2DM 05, n=1303" | (plot_id == "T2DM 05, n=1303" & all(fit_est <= 4)),
    plot_id != "T2DM 11, n=899" | (plot_id == "T2DM 11, n=899" & all(fit_est <= 6)),
    plot_id != "T2DM 16, n=791" | (plot_id == "T2DM 16, n=791" & all(fit_est <= 5))
  ) %>%
  ungroup() %>%
  ggplot(aes(x = time, group = iter)) +
  geom_step(
    data = ~ filter(.x, !is.na(estimate)), 
    aes(y = estimate), 
    linewidth = 0.5, 
    colour = "black"
  ) +
  geom_line(aes(y = fit_est, colour = dist), linewidth = 0.5, alpha = 0.02) +
  scale_colour_manual(
    values = c(
      "exp" = "blue",
      "gompertz" = "orange",
      "lnorm" = "purple",
      "llogis" = "deeppink",
      "weibull" = "red"
    )
  ) +
  scale_x_continuous(n.breaks = 4) +
  scale_y_continuous(n.breaks = 4) +
  labs(x = "Time in days", y = "Cumulative incidence of attrition") +
  facet_wrap(~plot_id, scales = "free") +
  theme_bw() +
  theme(
    text = element_text(size = 5),
    legend.position = "none"
  )

ggsave(
  "categorical_attrition/plots/suppl_l2f.png",
  plot_cind,
  dpi = 300,
  width = 2244/300,
  height = 1683/300,
  units = "in"
)

# Other/Miscellaneous (40 trials)
plot_cind <- cind_plot_df %>% 
  filter(cause == "Other/Miscellaneous") %>% 
  mutate(across(c(estimate, fit_est), ~ round((1 - .) * 100, 1))) %>% 
  group_by(plot_id, iter) %>%
  filter(
    plot_id != "COPD 10, n=906" | (plot_id == "COPD 10, n=906" & all(fit_est <= 5)),
    plot_id != "COPD 13, n=644" | (plot_id == "COPD 13, n=644" & all(fit_est <= 5)),
    plot_id != "T2DM 07, n=1162" | (plot_id == "T2DM 07, n=1162" & all(fit_est <= 2.5)),
    plot_id != "T2DM 14, n=807" | (plot_id == "T2DM 14, n=807" & all(fit_est <= 4)),
    plot_id != "T2DM 35, n=272" | (plot_id == "T2DM 35, n=272" & all(fit_est <= 10))
  ) %>%
  ungroup() %>%
  ggplot(aes(x = time, group = iter)) +
  geom_step(
    data = ~ filter(.x, !is.na(estimate)), 
    aes(y = estimate), 
    linewidth = 0.5, 
    colour = "black"
  ) +
  geom_line(aes(y = fit_est, colour = dist), linewidth = 0.5, alpha = 0.02) +
  scale_colour_manual(
    values = c(
      "exp" = "blue",
      "gompertz" = "orange",
      "lnorm" = "purple",
      "llogis" = "deeppink",
      "weibull" = "red"
    )
  ) +
  scale_x_continuous(n.breaks = 4) +
  scale_y_continuous(n.breaks = 4) +
  labs(x = "Time in days", y = "Cumulative incidence of attrition") +
  facet_wrap(~plot_id, scales = "free") +
  theme_bw() +
  theme(
    text = element_text(size = 5),
    legend.position = "none"
  )

ggsave(
  "categorical_attrition/plots/suppl_other.png",
  plot_cind,
  dpi = 300,
  width = 2244/300,
  height = 1683/300,
  units = "in"
)

# PI/Sponsor decision (6 trials)
plot_cind <- cind_plot_df %>% 
  filter(cause == "PI/Sponsor Decision") %>% 
  mutate(across(c(estimate, fit_est), ~ round((1 - .) * 100, 1))) %>% 
  group_by(plot_id, iter) %>%
  filter(
    plot_id != "BPH 01, n=1056" | (plot_id == "BPH 01, n=1056" & all(fit_est <= 2.5))
  ) %>%
  ungroup() %>%
  ggplot(aes(x = time, group = iter)) +
  geom_step(
    data = ~ filter(.x, !is.na(estimate)), 
    aes(y = estimate), 
    linewidth = 0.5, 
    colour = "black"
  ) +
  geom_line(aes(y = fit_est, colour = dist), linewidth = 0.5, alpha = 0.02) +
  scale_colour_manual(
    values = c(
      "exp" = "blue",
      "gompertz" = "orange",
      "lnorm" = "purple",
      "llogis" = "deeppink",
      "weibull" = "red"
    )
  ) +
  scale_x_continuous(n.breaks = 4) +
  scale_y_continuous(n.breaks = 4) +
  labs(x = "Time in days", y = "Cumulative incidence of attrition") +
  facet_wrap(~plot_id, scales = "free") +
  theme_bw() +
  theme(
    text = element_text(size = 5),
    legend.position = "none"
  )

ggsave(
  "categorical_attrition/plots/suppl_pi.png",
  plot_cind,
  dpi = 300,
  width = 2244/300,
  height = 1683/300,
  units = "in"
)

# Protocol violation (28 trials)
plot_cind <- cind_plot_df %>% 
  filter(cause == "Protocol Violation") %>% 
  mutate(across(c(estimate, fit_est), ~ round((1 - .) * 100, 1))) %>% 
  group_by(plot_id, iter) %>%
  filter(
    plot_id != "BPH 04, n=511" | (plot_id == "BPH 04, n=511" & all(fit_est <= 10)),
    plot_id != "Hyp 01, n=1039" | (plot_id == "Hyp 01, n=1039" & all(fit_est <= 5)),
    plot_id != "T2DM 24, n=566" | (plot_id == "T2DM 24, n=566" & all(fit_est <= 5))
  ) %>%
  ungroup() %>%
  ggplot(aes(x = time, group = iter)) +
  geom_step(
    data = ~ filter(.x, !is.na(estimate)), 
    aes(y = estimate), 
    linewidth = 0.5, 
    colour = "black"
  ) +
  geom_line(aes(y = fit_est, colour = dist), linewidth = 0.5, alpha = 0.02) +
  scale_colour_manual(
    values = c(
      "exp" = "blue",
      "gompertz" = "orange",
      "lnorm" = "purple",
      "llogis" = "deeppink",
      "weibull" = "red"
    )
  ) +
  scale_x_continuous(n.breaks = 4) +
  scale_y_continuous(n.breaks = 4) +
  labs(x = "Time in days", y = "Cumulative incidence of attrition") +
  facet_wrap(~plot_id, scales = "free") +
  theme_bw() +
  theme(
    text = element_text(size = 5),
    legend.position = "none"
  )

ggsave(
  "categorical_attrition/plots/suppl_pv.png",
  plot_cind,
  dpi = 300,
  width = 2244/300,
  height = 1683/300,
  units = "in"
)

# Voluntary withdrawal (64 trials)
plot_cind <- cind_plot_df %>% 
  filter(cause == "Voluntary Withdrawal") %>% 
  mutate(across(c(estimate, fit_est), ~ round((1 - .) * 100, 1))) %>% 
  group_by(plot_id, iter) %>%
  filter(
    plot_id != "Hyp 04, n=860" | (plot_id == "Hyp 04, n=860" & all(fit_est <= 5)),
    plot_id != "T2DM 04, n=1360" | (plot_id == "T2DM 04, n=1360" & all(fit_est <= 5)),
    plot_id != "T2DM 14, n=807" | (plot_id == "T2DM 14, n=807" & all(fit_est <= 4)),
    plot_id != "T2DM 35, n=272" | (plot_id == "T2DM 35, n=272" & all(fit_est <= 10))
  ) %>%
  ungroup() %>%
  ggplot(aes(x = time, group = iter)) +
  geom_step(
    data = ~ filter(.x, !is.na(estimate)), 
    aes(y = estimate), 
    linewidth = 0.5, 
    colour = "black"
  ) +
  geom_line(aes(y = fit_est, colour = dist), linewidth = 0.5, alpha = 0.02) +
  scale_colour_manual(
    values = c(
      "exp" = "blue",
      "gompertz" = "orange",
      "lnorm" = "purple",
      "llogis" = "deeppink",
      "weibull" = "red"
    )
  ) +
  scale_x_continuous(n.breaks = 4) +
  scale_y_continuous(n.breaks = 4) +
  labs(x = "Time in days", y = "Cumulative incidence of attrition") +
  facet_wrap(~plot_id, scales = "free") +
  theme_bw() +
  theme(
    text = element_text(size = 5),
    legend.position = "none"
  )

ggsave(
  "categorical_attrition/plots/suppl_vw.png",
  plot_cind,
  dpi = 300,
  width = 2244/300,
  height = 1683/300,
  units = "in"
)

