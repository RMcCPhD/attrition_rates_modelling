
# Load packages and functions
source("Survival_attrition/Scripts/00a_packages_ipd_functions.R")

# Create trial lookup ----------------------------------------------------------

# Import lookup from Extract_data
import_lkp_extracts <- read_csv("Extract_data/Data/lookup_trials_v2.csv")

# Summarise conditions (for harmonising below, run after creating)
sum_conditions <- create_lkp %>% 
  group_by(setting) %>% 
  reframe() %>% 
  mutate(
    condition = case_when(
      grepl("hypertension", setting) ~ "hypertension",
      grepl("pd", setting) & !grepl("copd", setting) ~ "parkinson's",
      grepl("benign", setting) ~ "benign prostatic hyperplasia",
      setting == "diabetic neuropathies" ~ "t2d",
      grepl("osterperosis", setting) ~ "osteoporosis",
      grepl("fibrosis", setting) ~ "pulmonary fibrosis",
      TRUE ~ setting
    )
  )

# Create attrition study trials lookup
create_lkp <- import_lkp_extracts %>% 
  filter(sponsor %in% c("boehringer", "eli lilly")) %>% 
  filter(include == 1) %>% 
  left_join(sum_conditions) %>% 
  select(ctgov:study, condition, duration, rand) %>% 
  mutate(
    across(
      trial:study,
      ~ gsub(",", "", .)
    )
  ) 

# Save
write_csv(
  create_lkp,
  "Survival_attrition/Data/lookup_trials.csv"
)
