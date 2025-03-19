
# Packages
library(tidyverse)
library(naniar)
library(readxl)
library(styler)
library(haven)
library(broom)
library(survival)
library(muhaz)
library(flexsurv)
library(gt)
library(RDP)
library(gridExtra)

# Prepare time-to-event data ---------------------------------------------------

get_data <- function() {
  
  # Assign dataset
  dat <- import_t2e %>% 
    select(ctgov:condition, usubjid:trtgrp, time:outcome) %>% 
    mutate(ctgov_outcome = paste0(ctgov, "_", outcome))
  
  # Combine ctgov and causes of attrition into a lookup
  ref_trials <- dat %>% 
    group_by(ctgov, outcome, ctgov_outcome) %>% 
    reframe() 
  
  # List individual trials by ctgov
  lst_trials_ref <- split(dat, dat$ctgov)
  
  # List individual trials by ctgov_outcome
  lst_trials <- split(dat, dat$ctgov_outcome)
  
  # Remove censored list datasets and those with <10 events
  lst_trials_removal <- mapply(function(df, name) {
    if (nrow(df) < 10 || str_detect(name, "Censor")) {
      return(NULL)
    }  else {
      return(df)
    }
  },
  df = lst_trials,
  name = names(lst_trials),
  SIMPLIFY = FALSE
  )
  
  lst_trials_removal <- lst_trials_removal[!sapply(lst_trials_removal, is.null)]
  
  # Get full trial dataset for each list item
  lst_trials_full <- lapply(names(lst_trials_removal), function(name) {
    replace_name <- names(lst_trials_ref)
    replace_name <- replace_name[sapply(replace_name, function(r_name) str_detect(name, r_name))]
    if(length(replace_name) > 0) {
      return(lst_trials_ref[[replace_name[1]]])
    } else {
      return(NULL)
    }
  })
  
  names(lst_trials_full) <- names(lst_trials_removal)
  
  # Identify events in each listed dataframe for the respective cause
  lst_trials_events <- lapply(seq_along(lst_trials_full), function(i) {
    lst_trials_full[[i]] %>%
      mutate(
        cause = gsub("^[^_]*_", "", names(lst_trials_full[i])),
        outcome_event = if_else(outcome %in% cause, 1, 0)
      ) %>% 
      select(-event) %>% 
      rename(event = outcome_event)
  })
  
  names(lst_trials_events) <- names(lst_trials_full)
  
  # Remove participants with t=0 (Censored participants who did not have an
  # end time in the IPD, could be completers or dropouts)
  for (i in seq_along(lst_trials_events)) {
    name <- names(lst_trials_events[i])
    lst_trials_events[[name]] <- lst_trials_events[[i]] %>% 
      filter(time != 0)
  }
  
  # Filter reference using final list names
  ref_trials <- ref_trials %>% 
    filter(ctgov_outcome %in% names(lst_trials_events))
  
  lst_output <- list(
    lst_coa_t2e = lst_trials_events,
    ref = ref_trials
  )
    
  return(lst_output)
}

# Fit parametric time-to-event models ------------------------------------------

fit_models <- function() {
  
  # Assign data
  dat <- get_data_outputs$lst_coa_t2e

  # Distributions used to fit models
  distr <- c("exp", "gengamma", "gompertz", "llogis", "lnorm", "weibull")

  # List to store survival models
  lst_flex <- vector("list", length = length(distr))

  # Function to fit models and handle errors
  fit_mdls <- function(distr, data) {
    tryCatch(
      {
        flexsurvreg(
          Surv(time, event) ~ 1, 
          dist = distr, 
          data = data
        )
      },
      error = function(e) {
        message(
          paste0(
            "Error fitting ", 
            paste0(unique(data$ctgov), "_", unique(data$cause)), 
            " ", 
            e
          )
        )
        return(NULL)
      }
    )
  }
  
  # Fit parametric models
  print("Fitting models")

  # Fit models for each distribution
  lst_flex <- lapply(distr, function(distr) {
    print(distr)
    lapply(dat, function(data) {
      fit_mdls(distr, data)
    })
  })
  
  # Remove NULL resulting from errors in fitting model
  lst_flex <- lapply(lst_flex, function(i) {
    Filter(Negate(is.null), i)
  })

  names(lst_flex) <- distr
  
  return(lst_flex)
  
}

# Get model outputs ------------------------------------------------------------

model_outputs <- function() {
  
  # Assign data
  dat <- get_fitted_models
  
  # Distributions used to fit models
  distr <- c("exp", "gengamma", "gompertz", "llogis", "lnorm", "weibull")
  
  print("Extracting flexsurv datasets")
  
  # Extract the datasets produced by flexsurv
  lst_model_data <- imap(
    get_fitted_models,
    ~ map(
      .x,
      ~ as_tibble(.x$data$Y) %>% select(time, status)
    )
  )
  
  print("Extracting model parameters")
  
  # Extract model parameters
  lst_params <- lapply(seq_along(dat), function(i) {
    lapply(seq_along(dat[[i]]), function(j) {
      t(dat[[i]][[j]][["res.t"]]) %>%
        as.data.frame() %>%
        rownames_to_column(var = "outputs") %>%
        filter(outputs %in% c("est", "se")) %>% 
        mutate(
          ctgov_cause = paste(names(dat[[i]][j])),
          dist = distr[i]
        ) %>% 
        separate(ctgov_cause, into = c("ctgov", "cause"), sep = "_")
    })
  })
  
  # Bind second-level lists
  lst_params <- lapply(seq_along(lst_params), function(i) {
    bind_rows(lst_params[[i]])
  })
  
  # Name top-level lists as respective distribution
  names(lst_params) <- distr
  
  print("Extracting performance metrics")
  
  # Extract and tidy performance metrics
  lst_fit <- lapply(seq_along(dat), function(i) {
    lapply(seq_along(dat[[i]]), function(j) {
      tibble(
        ctgov_cause = paste(names(dat[[i]][j])),
        dist = distr[i],
        loglik = dat[[i]][[j]][["loglik"]],
        aic = dat[[i]][[j]][["AIC"]]
      ) %>% 
        separate(ctgov_cause, into = c("ctgov", "cause"), sep = "_")
    })
  })
  
  # Bind second-level lists
  lst_fit <- lapply(seq_along(lst_fit), function(i) {
    bind_rows(lst_fit[[i]])
  })
  
  # Name top-level lists as respective distribution
  names(lst_fit) <- distr
  
  print("Extracting variance-covariance matrices")
  
  # Extract covariance matrices for model parameters
  lst_vcovs <- lapply(seq_along(dat), function(i) {
    lapply(seq_along(dat[[i]]), function(j) {
      dat[[i]][[j]][["cov"]]
    })
  })
  
  # Name second-level lists as respective trial_cause
  for (i in seq_along(lst_vcovs)) {
    names(lst_vcovs[[i]]) <- names(dat[[i]])
  }
  
  # Compile list of dataframes
  lst_output <- list(
    lst_t2e = lst_model_data,
    lst_params = lst_params,
    lst_fit = lst_fit,
    lst_cov = lst_vcovs
  )
  
  return(lst_output)
  
}
  
# Get model estimates ----------------------------------------------------------

model_estimates <- function() {
  
  # Assign data
  dat <- get_fitted_models
  dat2 <- get_model_outputs$lst_t2e
  
  # Distributions used to fit models
  distr <- c("exp", "gengamma", "gompertz", "llogis", "lnorm", "weibull")
  
  print("Getting model survival estimates")
  
  # Get parametric estimates
  lst_fitted_t2e <- imap(
    dat,
    ~ map(
      .x,
      ~summary(.x, type = "survival", ci = FALSE, tidy = TRUE)
    )
  )
  
  lst_fitted_t2e_tidy <- lapply(seq_along(lst_fitted_t2e), function(i) {
    lapply(seq_along(lst_fitted_t2e[[i]]), function(j) {
      lst_fitted_t2e[[i]][[j]] %>% 
        mutate(
          ctgov_cause = names(lst_fitted_t2e[[i]][j]),
          dist = names(lst_fitted_t2e[i])
        ) %>% 
        separate(
          ctgov_cause, 
          into = c("ctgov", "cause"), 
          sep = "_"
        )
    })
  })
  
  names(lst_fitted_t2e_tidy) <- distr
  
  for(i in seq_along(lst_fitted_t2e_tidy)) {
    names(lst_fitted_t2e_tidy[[i]]) <- names(dat[[i]])
  }
  
  print("Getting model hazard estimates")
  
  # Get hazard
  lst_fitted_hazards <- imap(
    dat,
    ~ map(
      .x,
      ~ summary(.x, type = "hazard", ci = FALSE, tidy = TRUE)
    )
  )
  
  lst_fitted_hazards_tidy <- lapply(seq_along(lst_fitted_hazards), function(i) {
    lapply(seq_along(lst_fitted_hazards[[i]]), function(j) {
      lst_fitted_hazards[[i]][[j]] %>% 
        mutate(
          ctgov_cause = names(lst_fitted_hazards[[i]][j]),
          dist = names(lst_fitted_hazards[i])
        ) %>% 
        separate(
          ctgov_cause, 
          into = c("ctgov", "cause"), 
          sep = "_"
        )
    })
  })
  
  names(lst_fitted_hazards_tidy) <- distr
  
  for(i in seq_along(lst_fitted_hazards_tidy)) {
    names(lst_fitted_hazards_tidy[[i]]) <- names(dat[[i]])
  }
  
  lst_haz_tidy <- lapply(seq_along(lst_fitted_hazards_tidy), function(i) {
    lapply(seq_along(lst_fitted_hazards_tidy[[i]]), function(j) {
      bind_rows(lst_fitted_hazards_tidy[[i]][[j]]) %>% 
        mutate(
          ctgov_cause = paste0(ctgov, "_", cause),
        ) %>% 
        select(-c(ctgov, cause)) %>% 
        arrange(time)
    })
  })
  
  names(lst_haz_tidy) <- distr
  
  for(i in seq_along(lst_haz_tidy)) {
    names(lst_haz_tidy[[i]]) <- names(dat[[i]])
  }
  
  # Compile list of dataframes
  lst_output <- list(
    lst_t2e = lst_fitted_t2e_tidy,
    lst_haz = lst_haz_tidy
  )
  
  return(lst_output)
  
}

# Get Kaplan-Meier estimates ---------------------------------------------------

km_estimates <- function() {
  
  # Assign data
  dat <- get_model_outputs$lst_t2e[[1]]
  dat2 <- get_model_estimates$lst_t2e
  dat3 <- get_model_estimates$lst_haz
  
  # Distributions used to fit models
  distr <- c("exp", "gengamma", "gompertz", "llogis", "lnorm", "weibull")
  
  # Function that gets an epsilon value that reduces data to 10-30 points
  get_epsilon <- function(data) {
    
    target_points <- c(10:28)
    epsilon <- 0.00005
    max_iterations = 1000
    step_size <- 0.00001
    
    for (i in 1:max_iterations) {
      
      rdp_data <- RamerDouglasPeucker(
        data$time,
        data$estimate,
        epsilon = epsilon
      )
      
      n_points <- nrow(rdp_data)
      
      if (n_points %in% target_points) {break} else  {
        epsilon <- epsilon + step_size
      }
    }
    
    return(epsilon)
    
  }
  
  print("Getting Kaplan-Meier estimates")
  
  # Get tidied Kaplan-Meier estimates
  lst_surv <- lapply(seq_along(dat), function(i) {
    tidy(
      survfit(
        Surv(dat[[i]]$time, dat[[i]]$status) ~ 1,
        data = dat[[i]]
      )
    )
  })
  
  # Reduce data points
  lst_surv_rdp <- lapply(seq_along(lst_surv), function(i) {
    epsilon <- get_epsilon(lst_surv[[i]])
    RamerDouglasPeucker(
      lst_surv[[i]]$time,
      lst_surv[[i]]$estimate,
      epsilon = epsilon
    )
  })
  
  names(lst_surv_rdp) <- names(dat)
  
  # Add descriptive variables
  lst_surv_rdp <- lapply(seq_along(lst_surv_rdp), function(i) {
    lst_surv_rdp[[i]] %>%
      rename(
        time = x, 
        estimate = y
      ) %>% 
      mutate(
        dist = "surv",
        ctgov_cause = names(lst_surv_rdp[i])
      ) %>% 
      separate(
        ctgov_cause, 
        into = c("ctgov", "cause"), 
        sep = "_"
      )
  })
  
  names(lst_surv_rdp) <- names(dat)
  
  # Bind rows to get tidy survival data for export
  df_surv <- bind_rows(lst_surv_rdp) %>% select(-dist)
  
  print("Joining KM to model estimates")
  
  # Add survival data into distribution lists
  lst_surv_rdp_listed <- vector("list", length = length(distr))
  names(lst_surv_rdp_listed) <- distr
  
  for(i in seq_along(lst_surv_rdp_listed)) {
    lst_surv_rdp_listed[[i]] <- lst_surv_rdp
  }
  
  # Join parametric estimates to survival and tidy
  lst_surv_tidy <- lapply(seq_along(dat2), function(i) {
    lapply(seq_along(dat2[[i]]), function(j) {
      bind_rows(
        dat2[[i]][[j]],
        lst_surv_rdp_listed[[i]][[j]]
      ) %>% 
        mutate(
          ctgov_cause = paste0(ctgov, "_", cause),
          estimate = if_else(!is.na(estimate), estimate, est)
        ) %>% 
        select(-c(est, ctgov, cause)) %>% 
        arrange(time) %>% 
        rename(type = dist) %>% 
        filter(
          time >= min(time[type == "surv"]) 
          & time <= max(time[type == "surv"])
        )
    })
  })
  
  names(lst_surv_tidy) <- distr
  
  for(i in seq_along(lst_surv_tidy)) {
    names(lst_surv_tidy[[i]]) <- names(dat2[[i]])
  }
  
  return(lst_surv_tidy)
  
}
  
# Generate time-to-event plots -------------------------------------------------

t2e_plots <- function() {
  
  # Assign data
  dat <- get_km_estimates
  
  # Distributions used to fit models
  distr <- c("exp", "gengamma", "gompertz", "llogis", "lnorm", "weibull")
  
  print("Generating plots")
  
  lst_plots <- lapply(seq_along(dat), function(i) {
    lapply(seq_along(dat[[i]]), function(j) {
      ggplot(
        data = dat[[i]][[j]],
        aes(
          x = time, 
          y = estimate
        )
      ) +
        geom_step(
          data = dat[[i]][[j]] %>% filter(type == "surv"),
          aes(
            x = time, 
            y = estimate
          ),
          linewidth = 2
        ) +
        geom_line(
          data = dat[[i]][[j]] %>% filter(type != "surv"),
          aes(
            x = time, 
            y = estimate
          ),
          linewidth = 2,
          colour = "salmon"
        )
    })
  })
  
  names(lst_plots) <- distr
  
  for(i in seq_along(lst_plots)) {
    names(lst_plots[[i]]) <- names(dat[[i]])
  }
  
  return(lst_plots)
  
}

# Generate hazard plots --------------------------------------------------------

haz_plots <- function() {
  
  # Assign data
  dat <- get_model_estimates$lst_haz
  
  # Distributions used to fit models
  distr <- c("exp", "gengamma", "gompertz", "llogis", "lnorm", "weibull")
  
  print("Generating plots")
  
  lst_plots <- lapply(seq_along(dat), function(i) {
    lapply(seq_along(dat[[i]]), function(j) {
      dat[[i]][[j]] %>% 
        ggplot(
          aes(
            x = time, 
            y = est
          )
        ) +
        geom_line(
          linewidth = 2,
          colour = "salmon"
        )
    })
  })
  
  names(lst_plots) <- distr
  
  for(i in seq_along(lst_plots)) {
    names(lst_plots[[i]]) <- names(dat[[i]])
  }
  
  return(lst_plots)
  
}
