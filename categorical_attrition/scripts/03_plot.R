
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


