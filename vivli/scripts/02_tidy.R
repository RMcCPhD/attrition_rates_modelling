
# Load packages and functions
source("Survival_attrition/Scripts/00a_packages_ipd_functions.R")

# Import lookup
import_lkp <- read_csv("Survival_attrition/Data/lookup_trials.csv") %>% 
  rename(studyid = study)

# Import IPD attrition extracts
import_boe <- readRDS("Survival_attrition/Data/boe.rds")
import_eli <- readRDS("Survival_attrition/Data/eli.rds")

## Prepare time-to-event dataset -----------------------------------------------

# Prepare boehringer trials
prep_boe <- import_lkp %>% 
  inner_join(import_boe) %>% 
  select(ctgov:sponsor, studyid:trtgrp, randt, disct, compl, reason) %>% 
  mutate(
    across(
      randt:disct,
      as.Date
    ),
    compl = as.integer(compl),
    time_days = as.numeric(disct - randt),
    time_week = time_days/7,
    event = case_when(
      time_week > duration | is.na(time_week) | time_week < 0 | compl == 1 ~ 0,
      TRUE ~ 1
    ),
    outcome = case_when(
      event == 0 ~ "Censored (Completer/Beyond intended duration/Time issue)",
      TRUE ~ reason
    ),
    time = case_when(
      event == 0 & (is.na(time_week) | time_week < 0) ~ 0,
      TRUE ~ time_days
    )
  ) %>% 
  select(ctgov:trtgrp, time, event, outcome, compl, reason) %>% 
  mutate(
    time = case_when(
      time == 0 & event != 0 ~ 1,
      time == 0 & event == 0 ~ 0,
      TRUE ~ time
    )
  )

# Prepare eli trials
prep_eli <- import_lkp %>% 
  inner_join(import_eli) %>% 
  select(ctgov:sponsor, studyid:trtgrp, randt, disct, compl, reason) %>% 
  mutate(
    across(
      randt:disct,
      as.Date
    ),
    compl = as.integer(compl),
    time_days = as.numeric(disct - randt),
    time_week = time_days/7,
    event = case_when(
      time_week > duration | is.na(time_week) | time_week < 0 | compl == 1 ~ 0,
      TRUE ~ 1
    ),
    time = case_when(
      event == 0 & (is.na(time_week) | time_week < 0) ~ 0,
      TRUE ~ time_days
    )
  ) %>% 
  select(ctgov:trtgrp, time, event, compl, reason) %>% 
  mutate(
    time = case_when(
      time == 0 & event != 0 ~ 1,
      time == 0 & event == 0 ~ 0,
      TRUE ~ time
    )
  )

# Bind together
prep_bind <- prep_boe %>% bind_rows(prep_eli)

# Harmonise causes of attrition
# Assume those missing completion status did not complete
# They are censored regardless because of missing times
t2e_lkp_causes <- prep_bind %>% 
  group_by(reason) %>% 
  reframe() %>% 
  mutate(
    cause = case_when(
      grepl("ADVERSE|Adv.|Adverse", reason) ~ "Adverse Event",
      grepl("LACK OF|Lack of", reason) ~ "Lack of Efficacy",
      grepl("LOST|Lost|lost", reason) ~ "Lost to Follow-up",
      grepl("Investigator|PHYSICIAN|Physician|Sponsor", reason) ~ "PI/Sponsor Decision",
      grepl("PROTOCOL|Protocol D|Protocol V", reason) ~ "Protocol Violation",
      grepl("^Patient|^Subject|^Voluntary|WITHDRAWAL|Withdrawal", reason) ~ "Voluntary Withdrawal",
      grepl("COMPLETED|Completed", reason) ~ "Completed",
      TRUE ~ "Other/Miscellaneous"
    )
  )

# Prepare time-to-event data
t2e_data <- prep_bind %>% 
  left_join(t2e_lkp_causes) %>% 
  mutate(
    outcome = case_when(
      event == 0 ~ "Censored (Completer/Beyond intended duration/Time issue)",
      TRUE ~ cause
    )
  ) %>% 
  select(ctgov:trtgrp, compl, cause, time, event, outcome) %>% 
  mutate(
    compl = if_else(is.na(compl), 0, compl)
  )

# Save
write_rds(
  t2e_data,
  "Survival_attrition/Data/t2e_all.rds"
)

## Check events and select trials ----------------------------------------------

# Check events per trial (where uncensored)
t2e_check_events <- t2e_data %>% 
  filter(!grepl("Censored", outcome)) %>% 
  group_by(ctgov, cause) %>%
  summarise(
    count = sum(event)
  ) %>%
  group_by(ctgov) %>%
  mutate(
    mt10 = sum(count >= 10)
  ) %>%
  pivot_wider(
    names_from = cause, 
    values_from = count,
    values_fill = 0
  )

# Select trials where at least 10 events occurred for at least one cause
# N=92 trials
t2e_select_trials <- t2e_check_events %>% 
  filter(mt10 != 0) %>% 
  select(ctgov)

# Trials included in analysis (N=92)
t2e_selected_trials <- t2e_data %>% filter(ctgov %in% t2e_select_trials$ctgov)

# Save
write_rds(
  t2e_selected_trials,
  "Survival_attrition/Data/t2e.rds"
)











