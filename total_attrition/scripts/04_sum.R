
# Produce summary statistics and tables
source("config.R")
library(tidyverse)

# Import data
imp_flow <- read_csv("vivli/data/summaries/participant_flow_formatted.csv")
imp_cind <- read_rds("total_attrition/processed_data/cind_sampled_df.rds")

# Factorise and order conditions, get values
flow_fmt <- imp_flow %>% 
  mutate(
    condition = factor(condition, levels = order_conditions),
    across(completers:`Voluntary Withdrawal`, ~ gsub("\\s*\\([^)]*\\)", "", .)),
    across(completers:dropouts,as.numeric)
  ) %>% 
  arrange(condition)

# Get person-years of follow-up and follow-up duration
py_fu <- imp_cind %>% 
  select(ctgov, time) %>% 
  distinct() %>% 
  group_by(ctgov) %>% 
  reframe(
    fu = max(time),
    fu_wk = round(fu/7, 1),
    fu_mth = round(fu/30, 1),
    fu_yr = round(fu/365, 1)
  )

# Supplementary table 2. Study characteristics #################################

# Trials, participants, follow-up, person-years of follow-up, events
tbl1 <- flow_fmt %>% 
  rename(
    n = participants,
    total_events = dropouts,
    ae_events = `Adverse Event`,
    loe_events = `Lack of Efficacy`,
    l2f_events = `Lost to Follow-up`,
    o_events = `Other/Miscellaneous`,
    pi_events = `PI/Sponsor Decision`,
    pv_events = `Protocol Violation`,
    vw_events = `Voluntary Withdrawal`
  ) %>% 
  left_join(py_fu %>% select(ctgov, fu_yr, fu_mth)) %>% 
  mutate(
    pyr = n * fu_yr,
    pyr_events = total_events/pyr
  ) %>% 
  select(ctgov:n, fu_mth, pyr, total_events:vw_events) %>% 
  arrange(condition, desc(n)) %>% 
  mutate(across(total_events:vw_events, ~ as.numeric(.)))

# Aggregated summary statistics for in-text results
tbl1 %>% 
  reframe(
    across(
      .cols = n:vw_events,
      .fns = list(
        median = ~ median(.x, na.rm = TRUE),
        min = ~ min(.x, na.rm = TRUE),
        max = ~ max(., na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    )
  )

# Save as publicly usable data
# write_csv(tbl1, "usable_data/total_attrition/study_chars.csv")

# Aggregate and format
tbl1_fmt <- tbl1 %>% 
  group_by(condition) %>% 
  reframe(
    n_trials = n(),
    across(
      .cols = n:vw_events,
      .fns = list(
        median = ~ median(.x, na.rm = TRUE),
        min = ~ min(.x, na.rm = TRUE),
        max = ~ max(., na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    )) %>% 
  mutate(
    n_trials = paste0(n_trials, " (", round(n_trials/90 * 100, 1), ")"),
    n = paste0(n_median, " [", n_min, ", ", n_max, "]"),
    fu_mth = paste0(fu_mth_median, " [", fu_mth_min, ", ", fu_mth_max, "]"),
    pyr = paste0(pyr_median, " [", pyr_min, ", ", pyr_max, "]"),
    total_events = paste0(total_events_median, " [", total_events_min, ", ", total_events_max, "]"),
    ae_events = paste0(ae_events_median, " [", ae_events_min, ", ", ae_events_max, "]"),
    loe_events = paste0(loe_events_median, " [", loe_events_min, ", ", loe_events_max, "]"),
    l2f_events = paste0(l2f_events_median, " [", l2f_events_min, ", ", l2f_events_max, "]"),
    o_events = paste0(o_events_median, " [", o_events_min, ", ", o_events_max, "]"),
    pi_events = paste0(pi_events_median, " [", pi_events_min, ", ", pi_events_max, "]"),
    pv_events = paste0(pv_events_median, " [", pv_events_min, ", ", pv_events_max, "]"),
    vw_events = paste0(vw_events_median, " [", vw_events_min, ", ", vw_events_max, "]"),
  ) %>% 
  select(condition, n_trials, n:vw_events)

# write_csv(tbl1_fmt, "total_attrition/summaries/tbl1.csv")