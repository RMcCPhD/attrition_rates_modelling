
rm(list = ls())
gc()

options(scipen = 999)
set.seed(123)

# Ordering of conditions for summary statistics
order_conditions <- c(
  "Type 2 diabetes", "COPD", "Hypertension", "Benign prostatic hyperplasia", 
  "Parkinson's disease", "Pulmonary fibrosis", "Asthma", "Osteoarthritis", 
  "Psoriasis", "Restless leg syndrome", "Osteoporosis"
)