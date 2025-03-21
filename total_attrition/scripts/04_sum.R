
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

# Table 1. Study characteristics
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
  select(ctgov:n, fu_mth, pyr, total_events:vw_events)

write_csv(tbl1, "total_attrition/summaries/tbl1.csv")
write_csv(tbl1, "usable_data/total_attrition/study_chars.csv")