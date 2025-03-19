
# Total attrition --------------------------------------------------------------

# Load functions
source("Survival_attrition/Scripts/00b_packages_ta_functions.R")

# Get names of saved outputs
z_filenames <- list.files("Survival_attrition/Output/Total/")

# Paste filepaths together
z_filepaths <- c()
for(i in seq_along(z_filenames)) {
  z_filepaths[i] <- paste0("Survival_attrition/Output/Total/", z_filenames[i])
}

# Import outputs
imports_lst <- list()
for (i in seq_along(z_filepaths)) {
  name <- gsub(".rds", "", z_filenames[i])
  imports_lst[[name]] <-   read_rds(z_filepaths[[i]])
}

# Tidy time-to-event data (semi-colon separated, row-collapsed values)
tidy_t2e <- bind_rows(imports_lst$get_km_estimates[[1]]) %>% 
  filter(type == "surv") %>%
  mutate(
    across(
      c(time, estimate),
      ~ formatC(.x, format = "e", digits = 3)
    )
  ) %>%
  select(
    ctgov,  everything()
  ) %>% 
  group_by(ctgov) %>% 
  summarise(
    dat = paste0(" t=", time, ", est=", estimate, " ", collapse = ";")
  ) %>% 
  ungroup()

# Get wide-format parameters and bind each distribution
data_parameters <- imports_lst$get_model_outputs$lst_params

tidy_parameters <- lapply(seq_along(data_parameters), function(i) {
  if (grepl("exp", names(data_parameters[i]))) {
    data_parameters[[i]] %>%
      pivot_longer(rate, names_to = "parameter") %>%
      pivot_wider(names_from = "outputs", values_from = "value")
  } else if (grepl("weibull|llogis", names(data_parameters[i]))) {
    data_parameters[[i]] %>%
      pivot_longer(c(shape, scale), names_to = "parameter") %>%
      pivot_wider(names_from = "outputs", values_from = "value")
  } else if (grepl("gengamma", names(data_parameters[i]))) {
    data_parameters[[i]] %>%
      pivot_longer(c(mu, sigma, Q), names_to = "parameter") %>%
      pivot_wider(names_from = "outputs", values_from = "value")
  } else if (grepl("lnorm", names(data_parameters[i]))) {
    data_parameters[[i]] %>%
      pivot_longer(c(meanlog, sdlog), names_to = "parameter") %>%
      pivot_wider(names_from = "outputs", values_from = "value")
  } else if (grepl("gompertz", names(data_parameters[i]))) {
    data_parameters[[i]] %>%
      pivot_longer(c(shape, rate), names_to = "parameter") %>%
      pivot_wider(names_from = "outputs", values_from = "value")
  }
})

tidy_parameters <- bind_rows(tidy_parameters) %>%
  mutate(
    across(
      c(est, se),
      ~ formatC(.x, format = "e", digits = 3)
    )
  )

# Bind rows of model fit
data_fit <- imports_lst$get_model_outputs$lst_fit

tidy_fit <- bind_rows(data_fit) %>%
  mutate(
    across(
      c(loglik, aic),
      ~ formatC(.x, format = "e", digits = 3)
    )
  )

# Convert vcov matrices to dataframes and tidy
data_covariance <- imports_lst$get_model_outputs$lst_cov
names(data_covariance) <- names(data_parameters)

tidy_covariance <- lapply(seq_along(data_covariance), function(i) {
  lapply(seq_along(data_covariance[[i]]), function(j) {
    as.data.frame(data_covariance[[i]][[j]]) %>%
      rownames_to_column(var = "parameter1") %>%
      mutate(
        ctgov = names(data_covariance[[i]][j]),
        dist = names(data_covariance[i])
      )
  })
})

names(tidy_covariance) <- names(data_covariance)

# Get wide-format data for each distribution with manual tibbles for any NULL
tidy_covariance <- lapply(seq_along(tidy_covariance), function(i) {
  lapply(seq_along(tidy_covariance[[i]]), function(j) {
    if (grepl("exp", names(tidy_covariance[i]))) {
      tidy_covariance[[i]][[j]] %>%
        rename(rate = V1) %>% 
        pivot_longer(
          rate,
          names_to = "parameter2",
          values_to = "r"
        )
      
    } else if (grepl("weibull|llogis", names(tidy_covariance[i]))) {
      tidy_covariance[[i]][[j]] %>%
        pivot_longer(
          c(shape, scale), 
          names_to = "parameter2", 
          values_to = "r"
        )
      
    } else if (
      grepl("gengamma", names(tidy_covariance[i])) 
      & sum(is.na(tidy_covariance[[i]][[j]])) == 0
    ) {
      tidy_covariance[[i]][[j]] %>%
        pivot_longer(
          c(mu, sigma, Q), 
          names_to = "parameter2", 
          values_to = "r"
        )
      
    } else if (grepl("lnorm", names(tidy_covariance[i]))) {
      tidy_covariance[[i]][[j]] %>%
        pivot_longer(
          c(meanlog, sdlog), 
          names_to = "parameter2", 
          values_to = "r"
        )
      
    } else if (grepl("gompertz", names(tidy_covariance[i]))) {
      tidy_covariance[[i]][[j]] %>%
        pivot_longer(
          c(shape, rate), 
          names_to = "parameter2", 
          values_to = "r"
        )
    }
  })
})

# Bind rows and select variables
tidy_covariance <- bind_rows(tidy_covariance) %>%
  select(ctgov, dist, parameter1, parameter2, r) %>%
  mutate(r = formatC(r, format = "e", digits = 3))

# Save outputs for export
lst_outputs_names <- grep("^tidy", ls(), value = TRUE, perl = TRUE)
lst_outputs_get <- lapply(lst_outputs_names, get)
names(lst_outputs_get) <- lst_outputs_names

# Save selected outputs
for (i in seq_along(lst_outputs_get)) {
  name <- names(lst_outputs_get[i])
  name <- gsub("tidy_", "", name)
  msg <- paste0("Saving ", name)
  print(msg)
  write_csv(
    lst_outputs_get[[i]],
    paste0("Survival_attrition/Export/Total/Model_outputs/", name, ".csv")
  )
}


# Causes of attrition ----------------------------------------------------------

# Load functions
source("Survival_attrition/Scripts/00c_packages_coa_functions.R")

# Get names of saved outputs
z_filenames <- list.files("Survival_attrition/Output/Causes/")

# Paste filepaths together
z_filepaths <- c()
for(i in seq_along(z_filenames)) {
  z_filepaths[i] <- paste0("Survival_attrition/Output/Causes/", z_filenames[i])
}

# Import outputs
imports_lst <- list()
for (i in seq_along(z_filepaths)) {
  name <- gsub(".rds", "", z_filenames[i])
  imports_lst[[name]] <-   read_rds(z_filepaths[[i]])
}

# Tidy time-to-event data (semi-colon separated, row-collapsed values)
tidy_t2e <- bind_rows(imports_lst$get_km_estimates[[1]]) %>% 
  filter(type == "surv") %>%
  mutate(
    across(
      c(time, estimate),
      ~ formatC(.x, format = "e", digits = 3)
    )
  ) %>%
  separate(
    ctgov_cause,
    into = c("ctgov", "cause"),
    sep = "_"
  ) %>% 
  select(
    ctgov, cause, everything()
  ) %>% 
  group_by(ctgov, cause) %>% 
  summarise(
    dat = paste0(" t=", time, ", est=", estimate, " ", collapse = ";")
  ) %>% 
  ungroup()

# Get wide-format parameters and bind each distribution
data_parameters <- imports_lst$get_model_outputs$lst_params

tidy_parameters <- lapply(seq_along(data_parameters), function(i) {
  if (grepl("exp", names(data_parameters[i]))) {
    data_parameters[[i]] %>%
      pivot_longer(rate, names_to = "parameter") %>%
      pivot_wider(names_from = "outputs", values_from = "value")
  } else if (grepl("weibull|llogis", names(data_parameters[i]))) {
    data_parameters[[i]] %>%
      pivot_longer(c(shape, scale), names_to = "parameter") %>%
      pivot_wider(names_from = "outputs", values_from = "value")
  } else if (grepl("gengamma", names(data_parameters[i]))) {
    data_parameters[[i]] %>%
      pivot_longer(c(mu, sigma, Q), names_to = "parameter") %>%
      pivot_wider(names_from = "outputs", values_from = "value")
  } else if (grepl("lnorm", names(data_parameters[i]))) {
    data_parameters[[i]] %>%
      pivot_longer(c(meanlog, sdlog), names_to = "parameter") %>%
      pivot_wider(names_from = "outputs", values_from = "value")
  } else if (grepl("gompertz", names(data_parameters[i]))) {
    data_parameters[[i]] %>%
      pivot_longer(c(shape, rate), names_to = "parameter") %>%
      pivot_wider(names_from = "outputs", values_from = "value")
  }
})

tidy_parameters <- bind_rows(tidy_parameters) %>%
  mutate(
    across(
      c(est, se),
      ~ formatC(.x, format = "e", digits = 3)
    )
  )

# Bind rows of model fit
data_fit <- imports_lst$get_model_outputs$lst_fit

tidy_fit <- bind_rows(data_fit) %>%
  mutate(
    across(
      c(loglik, aic),
      ~ formatC(.x, format = "e", digits = 3)
    )
  )

# Convert vcov matrices to dataframes and tidy
data_covariance <- imports_lst$get_model_outputs$lst_cov
names(data_covariance) <- names(data_parameters)

tidy_covariance <- lapply(seq_along(data_covariance), function(i) {
  lapply(seq_along(data_covariance[[i]]), function(j) {
    as.data.frame(data_covariance[[i]][[j]]) %>%
      rownames_to_column(var = "parameter1") %>%
      mutate(
        ctgov_id = names(data_covariance[[i]][j]),
        dist = names(data_covariance[i])
      ) %>%
      separate(ctgov_id, into = c("ctgov", "cause"), sep = "_")
  })
})

names(tidy_covariance) <- names(data_covariance)

# Get wide-format data for each distribution with manual tibbles for any NULL
tidy_covariance <- lapply(seq_along(tidy_covariance), function(i) {
  lapply(seq_along(tidy_covariance[[i]]), function(j) {
    if (grepl("exp", names(tidy_covariance[i]))) {
      tidy_covariance[[i]][[j]] %>%
        rename(rate = V1) %>% 
        pivot_longer(
          rate,
          names_to = "parameter2",
          values_to = "r"
        )
      
    } else if (grepl("weibull|llogis", names(tidy_covariance[i]))) {
      tidy_covariance[[i]][[j]] %>%
        pivot_longer(
          c(shape, scale), 
          names_to = "parameter2", 
          values_to = "r"
        )
      
    } else if (
      grepl("gengamma", names(tidy_covariance[i])) 
      & sum(is.na(tidy_covariance[[i]][[j]])) == 0
    ) {
      tidy_covariance[[i]][[j]] %>%
        pivot_longer(
          c(mu, sigma, Q), 
          names_to = "parameter2", 
          values_to = "r"
        )
      
    } else if (grepl("lnorm", names(tidy_covariance[i]))) {
      tidy_covariance[[i]][[j]] %>%
        pivot_longer(
          c(meanlog, sdlog), 
          names_to = "parameter2", 
          values_to = "r"
        )
      
    } else if (grepl("gompertz", names(tidy_covariance[i]))) {
      tidy_covariance[[i]][[j]] %>%
        pivot_longer(
          c(shape, rate), 
          names_to = "parameter2", 
          values_to = "r"
        )
    }
  })
})

# Bind rows and select variables
tidy_covariance <- bind_rows(tidy_covariance) %>%
  select(ctgov, cause, dist, parameter1, parameter2, r) %>%
  mutate(r = formatC(r, format = "e", digits = 3))

# Save outputs for export
lst_outputs_names <- grep("^tidy", ls(), value = TRUE, perl = TRUE)
lst_outputs_get <- lapply(lst_outputs_names, get)
names(lst_outputs_get) <- lst_outputs_names

# Save selected outputs
for (i in seq_along(lst_outputs_get)) {
  name <- names(lst_outputs_get[i])
  name <- gsub("tidy_", "", name)
  msg <- paste0("Saving ", name)
  print(msg)
  write_csv(
    lst_outputs_get[[i]],
    paste0("Survival_attrition/Export/Causes/Model_outputs/", name, ".csv")
  )
}

