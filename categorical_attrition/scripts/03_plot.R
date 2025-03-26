
# Plot estimates from empirical data and models
# Plot files are not saved in repo to save file space, but can be obtained
# by running through the scripts

source("config.R")
source("categorical_attrition/scripts/00_functions.R")
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
  width = 190,
  height = 160,
  units = "mm"
)

### Compare best and second best ###############################################

# Get model metrics for best and second-best
best_2ndbest <- read_csv("categorical_attrition/data/fit.csv") %>% 
  filter(dist != "gengamma" & !(ctgov %in% c("NCT01131676", "NCT00274573"))) %>% 
  group_by(ctgov, cause) %>% 
  slice_min(aic, n = 2) %>% 
  slice_head(n = 2) |> 
  ungroup() %>% 
  left_join(imp_params) %>% 
  select(-loglik, -aic, -se) %>% 
  pivot_wider(names_from = "parameter", values_from = "est")

ci_cind <- lst_prep$cind %>% 
  unnest(cols = "cind") %>% 
  group_by(ctgov, cause) %>% 
  mutate(time_seq = list(seq(min(time), max(time), by = 1))) %>% 
  unnest(time_seq) %>% 
  select(ctgov, cause, time_seq) %>% 
  left_join(
    lst_prep$cind %>% 
      unnest(cols = "cind") %>% 
      rename(time_seq = time)
  ) %>% 
  distinct()

best_2ndbest_mdls <- ci_cind %>% 
  filter(ctgov %in% best_2ndbest$ctgov & cause %in% best_2ndbest$cause) %>% 
  left_join(best_2ndbest, by = c("ctgov", "cause")) %>% 
  group_by(ctgov, cause) %>%
  mutate(
    mdl = pmap_dbl(
      list(time_seq, dist, rate, shape, scale, meanlog, sdlog),
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
    )
  )

mdl_df <- ind_cond_n %>% 
  left_join(
    best_2ndbest_mdls %>% 
      select(ctgov, cause, time_seq, estimate, dist, mdl) 
  )

# Adverse event
supp_ae <- mdl_df %>% 
  filter(cause == "Adverse Event") %>% 
  ggplot(aes(x = time_seq)) +
  geom_step(
    data = ~ filter(.x, !is.na(estimate)),
    aes(y = 1 - estimate),
    linewidth = 0.5,
    colour = "black"
  ) +
  geom_line(aes(y = 1 - mdl, colour = dist), linewidth = 0.5) +
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
    legend.position = "none"
  )

ggsave(
  "categorical_attrition/plots/suppl_ae.png",
  supp_ae,
  dpi = 300,
  width = 4488/300,
  height = 2524/300,
  units = "in"
)

# Lack of Efficacy
supp_loe <- mdl_df %>% 
  filter(cause == "Lack of Efficacy") %>% 
  ggplot(aes(x = time_seq)) +
  geom_step(
    data = ~ filter(.x, !is.na(estimate)),
    aes(y = 1 - estimate),
    linewidth = 0.5,
    colour = "black"
  ) +
  geom_line(aes(y = 1 - mdl, colour = dist), linewidth = 0.5) +
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
    legend.position = "none"
  )

ggsave(
  "categorical_attrition/plots/suppl_loe.png",
  supp_loe,
  dpi = 300,
  width = 4488/300,
  height = 2524/300,
  units = "in"
)

# Lost to Follow-up
supp_l2f <- mdl_df %>% 
  filter(cause == "Lost to Follow-up") %>% 
  ggplot(aes(x = time_seq)) +
  geom_step(
    data = ~ filter(.x, !is.na(estimate)),
    aes(y = 1 - estimate),
    linewidth = 0.5,
    colour = "black"
  ) +
  geom_line(aes(y = 1 - mdl, colour = dist), linewidth = 0.5) +
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
    legend.position = "none"
  )

ggsave(
  "categorical_attrition/plots/suppl_l2f.png",
  supp_l2f,
  dpi = 300,
  width = 4488/300,
  height = 2524/300,
  units = "in"
)

# Other/Miscellaneous
supp_o <- mdl_df %>% 
  filter(cause == "Other/Miscellaneous") %>% 
  ggplot(aes(x = time_seq)) +
  geom_step(
    data = ~ filter(.x, !is.na(estimate)),
    aes(y = 1 - estimate),
    linewidth = 0.5,
    colour = "black"
  ) +
  geom_line(aes(y = 1 - mdl, colour = dist), linewidth = 0.5) +
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
    legend.position = "none"
  )

ggsave(
  "categorical_attrition/plots/suppl_o.png",
  supp_o,
  dpi = 300,
  width = 4488/300,
  height = 2524/300,
  units = "in"
)

# PI/Sponsor Decision
supp_pi <- mdl_df %>% 
  filter(cause == "PI/Sponsor Decision") %>% 
  ggplot(aes(x = time_seq)) +
  geom_step(
    data = ~ filter(.x, !is.na(estimate)),
    aes(y = 1 - estimate),
    linewidth = 0.5,
    colour = "black"
  ) +
  geom_line(aes(y = 1 - mdl, colour = dist), linewidth = 0.5) +
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
    legend.position = "none"
  )

ggsave(
  "categorical_attrition/plots/suppl_pi.png",
  supp_pi,
  dpi = 300,
  width = 4488/300,
  height = 2524/300,
  units = "in"
)

# Protocol Violation
supp_pv <- mdl_df %>% 
  filter(cause == "Protocol Violation") %>% 
  ggplot(aes(x = time_seq)) +
  geom_step(
    data = ~ filter(.x, !is.na(estimate)),
    aes(y = 1 - estimate),
    linewidth = 0.5,
    colour = "black"
  ) +
  geom_line(aes(y = 1 - mdl, colour = dist), linewidth = 0.5) +
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
    legend.position = "none"
  )

ggsave(
  "categorical_attrition/plots/suppl_pv.png",
  supp_pv,
  dpi = 300,
  width = 4488/300,
  height = 2524/300,
  units = "in"
)

# Voluntary Withdrawal
supp_vw <- mdl_df %>% 
  filter(cause == "Voluntary Withdrawal") %>% 
  ggplot(aes(x = time_seq)) +
  geom_step(
    data = ~ filter(.x, !is.na(estimate)),
    aes(y = 1 - estimate),
    linewidth = 0.5,
    colour = "black"
  ) +
  geom_line(aes(y = 1 - mdl, colour = dist), linewidth = 0.5) +
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
    legend.position = "none"
  )

ggsave(
  "categorical_attrition/plots/suppl_vw.png",
  supp_vw,
  dpi = 300,
  width = 4488/300,
  height = 2524/300,
  units = "in"
)
