
# Load packages and functions
source("Survival_attrition/Scripts/00a_packages_ipd_functions.R")

# Import time-to-event data for selected 92 trials
import_t2e <- readRDS("Survival_attrition/Data/t2e.rds")

# Participant flow -------------------------------------------------------------

# Participants (randomised, dropouts, completers)
sum_trial_participants <- import_t2e %>% 
  group_by(ctgov, condition) %>% 
  reframe(
    rand = length(ctgov),
    drop = sum(compl == 0),
    compl = sum(compl == 1)
  )

# Causes of attrition
sum_trial_causes <- import_t2e %>% 
  group_by(ctgov, cause) %>% 
  reframe(
    count = n()
  )

# Check that sum of causes is equal to participant values
# Both checks should be 1
sum_trial_check_flow <- sum_trial_participants %>% 
  left_join(sum_trial_causes) %>% 
  group_by(ctgov) %>% 
  mutate(
    check_count_drop = sum(count[cause != "Completed"]),
    check_count_compl = sum(count[cause == "Completed"])
  ) %>% 
  select(ctgov:compl, contains("check")) %>% 
  distinct() %>% 
  mutate(
    match_drop = if_else(drop == check_count_drop, 1, 0),
    match_compl = if_else(compl == check_count_compl, 1, 0)
  )
  
# Join causes to participants
sum_trial_flow <- sum_trial_participants %>% 
  left_join(sum_trial_causes) %>% 
  pivot_wider(
    names_from = cause,
    values_from = count,
    values_fill = 0
  ) %>% 
  select(-`Completed`)

# Participant flow formatted for tabulation
sum_trial_flow_fmt <- sum_trial_flow %>% 
  pivot_longer(
    `Adverse Event`:`PI/Sponsor Decision`,
    names_to = "cause",
    values_to = "count"
  ) %>% 
  mutate(
    drop_pcnt = round(drop/rand * 100, 2),
    compl_pcnt = round(compl/rand * 100, 2),
    completers = paste0(compl, " (", compl_pcnt, ")"),
    dropouts = paste0(drop, " (", drop_pcnt, ")"),
    pcnt = round(count/rand * 100, 2),
    n_pcnt = paste0(count, " (", pcnt, ")")
  ) %>% 
  select(
    ctgov, condition, participants = rand,
    completers, dropouts, cause, n_pcnt
  ) %>%
  pivot_wider(
    names_from = cause,
    values_from = n_pcnt,
    values_fill = "0 (0)"
  ) %>%
  select(ctgov:dropouts, everything())

# Save summarised trial-level participant flow
write_csv(
  sum_trial_flow,
  "Survival_attrition/Summaries/participant_flow_raw.csv"
)

write_csv(
  sum_trial_flow_fmt,
  "Survival_attrition/Summaries/participant_flow_formatted.csv"
)

# Aggregated summary statistics ------------------------------------------------

# Participants and total attrition
sum_aggregate <- import_t2e %>% 
  group_by(ctgov) %>% 
  reframe(
    participants = length(ctgov),
    drp_pcnt = (sum(compl == 0)/participants) * 100
  ) %>% 
  mutate(
    mean_attr = round(mean(drp_pcnt), 2),
    sd_attr = round(sd(drp_pcnt), 2),
    median_attr = round(median(drp_pcnt), 2),
    iqr_attr = round(IQR(drp_pcnt), 2),
    min_attr = round(min(drp_pcnt)),
    max_attr = round(max(drp_pcnt)),
    mean_n = round(mean(participants), 2),
    sd_n = round(sd(participants), 2),
    median_n = round(median(participants), 2),
    iqr_n = round(IQR(participants), 2),
    min_n = round(min(participants)),
    max_n = round(max(participants))
  ) %>% 
  select(-ctgov, -drp_pcnt, -participants) %>%
  distinct() %>%
  pivot_longer(
    mean_attr:max_n, 
    names_to = "statistic", 
    values_to = "value"
  )

# Causes of attrition
sum_aggregate_causes <- import_t2e %>% 
  group_by(cause) %>%
  reframe(count = n()) %>% 
  filter(cause != "Completed") %>%
  mutate(pcnt = round(count/sum(count) * 100, 2)) %>% 
  arrange(desc(pcnt))

# Conditions
sum_aggregate_conditions <- import_t2e %>% 
  group_by(ctgov, condition) %>%
  reframe(freq = n()) %>%
  mutate(freq = 1) %>%
  group_by(condition) %>%
  reframe(
    freq = sum(freq),
    pcnt = freq/92 * 100
  ) %>%
  ungroup() %>%
  arrange(desc(pcnt))

# Save aggregated summaries
write_csv(
  sum_aggregate,
  "Survival_attrition/Summaries/aggregated_participants_attrition.csv"
)

write_csv(
  sum_aggregate_causes,
  "Survival_attrition/Summaries/aggregated_attrition_causes.csv"
)

write_csv(
  sum_aggregate_conditions,
  "Survival_attrition/Summaries/aggregated_conditions.csv"
)










