
### Total attrition ------------------------------------------------------------

# Load functions
source("Survival_attrition/Scripts/00b_packages_ta_functions.R")

# Import time-to-event dataset
import_t2e <- readRDS("Survival_attrition/Data/t2e.rds")

# Prepare time-to-event data
get_data_outputs <- get_data()

# Fit parametric time-to-event models
get_fitted_models <- fit_models()

# Get model outputs
get_model_outputs <- model_outputs()

# Get model estimates
get_model_estimates <- model_estimates()

# Get Kaplan-Meier estimates
get_km_estimates <- km_estimates()

# Generate time-to-event plots
get_t2e_plots <- t2e_plots()

# Generate hazard plots
get_haz_plots <- haz_plots()

# Get names of output objects that are lists and not of plots
lst_outputs <- ls()[sapply(ls(), function(x) is.list(get(x)))]
lst_outputs_names <- grep("^get(?!.*models$)", lst_outputs, value = TRUE, perl = TRUE)
lst_outputs_get <- lapply(lst_outputs_names, get)
names(lst_outputs_get) <- lst_outputs_names

# Save selected outputs
for (i in seq_along(lst_outputs_get)) {
  name <- names(lst_outputs_get[i])
  msg <- paste0("Saving ", name)
  print(msg)
  write_rds(
    lst_outputs_get[[i]],
    paste0("Survival_attrition/Output/Total/", name, ".rds")
  )
}

### Causes of attrition --------------------------------------------------------

# Load functions
source("Survival_attrition/Scripts/00c_packages_coa_functions.R")

# Import time-to-event dataset
import_t2e <- readRDS("Survival_attrition/Data/t2e.rds")

get_data_outputs <- get_data()
get_fitted_models <- fit_models()
get_model_outputs <- model_outputs()
get_model_estimates <- model_estimates()
get_km_estimates <- km_estimates()

get_t2e_plots <- t2e_plots()
get_haz_plots <- haz_plots()

lst_outputs <- ls()[sapply(ls(), function(x) is.list(get(x)))]
lst_outputs_names <- grep("^get(?!.*models$)", lst_outputs, value = TRUE, perl = TRUE)
lst_outputs_get <- lapply(lst_outputs_names, get)
names(lst_outputs_get) <- lst_outputs_names

for (i in seq_along(lst_outputs_get)) {
  name <- names(lst_outputs_get[i])
  msg <- paste0("Saving ", name)
  print(msg)
  write_rds(
    lst_outputs_get[[i]],
    paste0("Survival_attrition/Output/Causes/", name, ".rds")
  )
}

