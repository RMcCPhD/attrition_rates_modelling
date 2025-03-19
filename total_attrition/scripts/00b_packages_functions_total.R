
# Packages
library(styler)
library(tidyverse)
library(flexsurv)
library(MASS)
library(zoo)
library(cowplot)
library(audio)
library(beepr)
library(gt)
library(RColorBrewer)

# Turn off scientific notation
options(scipen = 999)

# Get transformed datasets -----------------------------------------------------

get_data <- function() {
  
  dat <- import_t2e
  dat2 <- import_parameters
  dat3 <- import_fit
  dat4 <- import_covariance
  
  print("Transforming Kaplan-Meier data")
  
  # Transform Kaplan-Meier data
  t2e_transform <- dat %>% 
    group_by(ctgov) %>% 
    separate_rows(
      dat,
      sep = " ; "
    ) %>% 
    separate(
      dat,
      into = c("time", "estimate"),
      sep = ", "
    ) %>% 
    ungroup() %>% 
    mutate(
      across(
        c(time, estimate),
        ~ gsub("t|est|\\=", "", .)
      ),
      across(
        c(time, estimate),
        ~ as.numeric(.)
      )
    )
  
  # Nest Kaplan-Meier data
  t2e <- t2e_transform %>% 
    group_by(ctgov) %>% 
    nest(.key = "t2e") %>% 
    ungroup()
  
  # Join parameters and fit metrics
  params_fit <- dat2 %>% left_join(dat3) %>% arrange(ctgov)
  
  print("Getting best models")
  
  # Select the best 3 models per cause and trial based on AIC
  # Many cases where 4 are selected because the AIC values were the same between
  # two or more
  get_best_models <- params_fit %>% 
    dplyr::select(ctgov, dist, aic) %>% 
    filter(dist != "gengamma") %>% 
    distinct() %>% 
    group_by(ctgov) %>% 
    slice_min(aic, n = 1) %>% 
    slice_head(n = 1) %>% 
    ungroup()
  
  print("Nesting best model parameters")
  
  # Nest parameters for best performing models
  best_models <- params_fit %>% 
    inner_join(
      get_best_models %>% 
        dplyr::select(ctgov, dist)
    ) %>% 
    dplyr::select(-aic, -loglik) %>% 
    nest(
      .key = "nested_params",
      .by = c("ctgov")
    ) %>% 
    mutate(
      nested_params = 
        map(
          nested_params,
          ~ .x %>%
            dplyr::select(-se) %>% 
            group_by(dist) %>% 
            pivot_wider(names_from = parameter, values_from = est) %>% 
            nest(.key = "params") %>% 
            ungroup() %>% 
            arrange(dist) %>%
            mutate(
              params = 
                map(
                  params,
                  ~ as.matrix(
                    .x %>% 
                      select_if(~ !any(is.na(.))) %>% 
                      drop_na() %>% 
                      dplyr::select(order(colnames(.)))
                  )
                )
            )
        )
    ) %>% 
    arrange(ctgov)
  
  print("Nesting variance-covariance")
  
  # Nest variance-covariance matrices
  vcov <- import_covariance %>% 
    inner_join(
      get_best_models %>% 
        dplyr::select(ctgov, dist)
    ) %>% 
    nest(
      .key = "distr_cov",
      .by = c("ctgov")
    ) %>% 
    mutate(
      distr_cov = 
        map(
          distr_cov,
          ~ .x %>% 
            group_by(dist) %>% 
            spread(., key = parameter2, value = r) %>% 
            nest(.key = "vcov") %>% 
            ungroup() %>%  
            mutate(
              vcov = map(
                vcov, 
                ~ {
                  mat <- as.matrix(
                    .x %>% 
                      select_if(~ !any(is.na(.))) %>% 
                      drop_na() %>% 
                      dplyr::select(-parameter1) %>% 
                      dplyr::select(order(colnames(.)))
                  )
                  rownames(mat) <- colnames(mat)
                  mat
                }
              )
            )
        )
    ) %>% 
    arrange(ctgov)
  
  lst_output <- list(
    df_t2e = t2e_transform,
    nested_t2e = t2e,
    parameters_fit = params_fit,
    best_models = best_models,
    covariance = vcov
  )
  
  return(lst_output)
  
}

# Sample parameters from a multivariate normal distribution --------------------

get_samples <- function() {
  
  dat <- data$nested_t2e
  dat2 <- data$best_models
  dat3 <- data$covariance
  
  # List to store samples for each model
  lst_samples <- vector("list", length = nrow(dat))
  names(lst_samples) <- dat$ctgov
  
  for(i in seq_along(dat2$nested_params)) {
    
    name_loop <- dat$ctgov[i]
    
    for(j in seq_along(dat2$nested_params[[i]]$params)) {
      
      dist <- paste0(dat2$nested_params[[i]]$dist[[j]])
      surv <- dat$t2e[[i]]
      
      # mvrnorm arguments
      n = 100
      mu = dat2$nested_params[[i]]$params[[j]]  # Estimate
      sigma = dat3$distr_cov[[i]]$vcov[[j]]  # Covariance matrix
      
      # Sample parameters with mvrnorm
      tryCatch(
        {
          lst_samples[[i]][[dist]] <- mvrnorm(n = n, mu = mu, Sigma = sigma)
        },
        error = function(e) {
          message(paste0("Error sampling ", name_loop, " ", dist, " Error: ", e))
          return(NULL)
        }
      )
    }
  }
  
  return(lst_samples)
  
}

# Generate parametric estimates using sampled parameters -----------------------

get_estimates <- function() {
  
  dat <- data$nested_t2e
  dat2 <- sample_parameters
  
  # List to store estimates for each model
  lst_estimates <- vector("list", length = nrow(dat))
  names(lst_estimates) <- dat$ctgov
  
  for(i in seq_along(dat2)) {
    
    surv <- dat$t2e[[i]] # Non-parametric survival data
    name_loop <- names(dat2[i])  # ctgov_cause
    print(name_loop)  # Check for where errors occur
    
    for(j in seq_along(dat2[[i]])) {
      
      time <- round(seq(min(surv$time), max(surv$time)))
      dist <- names(dat2[[i]][j])  # Distribution name
      print(dist)  # Check for where errors occur
      
      for(k in seq_along(dat2[[i]][[j]])) {
        
        # Conditionally estimate fitted survival
        if(dist == "exp") {
          
          name_k <- paste0("iter_", k)  # Iteration
          rate <- dat2[[i]][[j]][[k]]  # Sampled parameter
          
          fit_est <- 1 - pexp(  # Fitted estimates
            q = time, 
            rate = exp(rate)  # Inverse log
          )
          
          # Store as data frame with relevant variables
          lst_estimates[[name_loop]][[dist]][[name_k]] <- data.frame(
            iter = k,
            fit_est = fit_est,
            time = time
          ) %>% 
            left_join(
              surv,
              by = "time"
            )
          
        } else if (dist == "gengamma") {
          
          params <- as.data.frame(dat2[[i]][[j]])
          
          for(n in 1:nrow(params)) {
            
            name_k <- paste0("iter_", n)
            
            mu_value <- params[n, "mu"]
            sigma_value <- params[n, "sigma"]  # Inverse log
            Q_value <- params[n, "Q"]
            
            fit_est <- 1 - pgengamma(
              q = time,
              mu = mu_value,
              sigma = exp(sigma_value),
              Q = Q_value
            )
            
            lst_estimates[[name_loop]][[dist]][[name_k]] <- data.frame(
              iter = n,
              fit_est = fit_est,
              time = time
            ) %>% 
              left_join(
                surv,
                by = "time"
              )
          }
          
        } else if (dist == "gompertz") {
          
          params <- as.data.frame(dat2[[i]][[j]])
          
          for(n in 1:nrow(params)) {
            
            name_k <- paste0("iter_", n)
            
            rate_value <- params[n, "rate"]
            shape_value <- params[n, "shape"]
            
            fit_est <- 1 - pgompertz(
              q = time,
              rate = exp(rate_value),  # Inverse log
              shape = shape_value
            )
            
            lst_estimates[[name_loop]][[dist]][[name_k]] <- data.frame(
              iter = n,
              fit_est = fit_est,
              time = time
            ) %>% 
              left_join(
                surv,
                by = "time"
              )
          }
          
        } else if (dist == "llogis") {
          
          params <- as.data.frame(dat2[[i]][[j]])
          
          for(n in 1:nrow(params)) {
            
            name_k <- paste0("iter_", n)
            
            scale_shape <- params[n, "scale"]
            shape_value <- params[n, "shape"]
            
            fit_est <- 1 - pllogis(
              q = time,
              scale = exp(scale_shape),  # Inverse log
              shape = exp(shape_value)  # Inverse log
            )
            
            lst_estimates[[name_loop]][[dist]][[name_k]] <- data.frame(
              iter = n,
              fit_est = fit_est,
              time = time
            ) %>% 
              left_join(
                surv,
                by = "time"
              )
          }
          
        } else if (dist == "lnorm") {
          
          params <- as.data.frame(dat2[[i]][[j]])
          
          for(n in 1:nrow(params)) {
            
            name_k <- paste0("iter_", n)
            
            meanlog_shape <- params[n, "meanlog"]
            sdlog_value <- params[n, "sdlog"]
            
            fit_est <- 1 - plnorm(
              q = time,
              meanlog = meanlog_shape,
              sdlog = exp(sdlog_value)  # Inverse log
            )
            
            lst_estimates[[name_loop]][[dist]][[name_k]] <- data.frame(
              iter = n,
              fit_est = fit_est,
              time = time
            ) %>% 
              left_join(
                surv,
                by = "time"
              )
          }
          
        } else if (dist == "weibull") {
          
          params <- as.data.frame(dat2[[i]][[j]])
          
          for(n in 1:nrow(params)) {
            
            name_k <- paste0("iter_", n)
            
            scale_shape <- params[n, "scale"]
            shape_value <- params[n, "shape"]
            
            fit_est <- 1 - pweibull(
              q = time,
              scale = exp(scale_shape),  # Inverse log
              shape = exp(shape_value)  # Inverse log
            )
            
            lst_estimates[[name_loop]][[dist]][[name_k]] <- data.frame(
              iter = n,
              fit_est = fit_est,
              time = time
            ) %>% 
              left_join(
                surv,
                by = "time"
              )
          }
        }
      }
    }
  }
  
  return(lst_estimates)
  
}

# Get hazard estimates from model parameters -----------------------------------

get_haz <- function() {
  
  # Conditionally import data
  dat <- data$nested_t2e
  dat2 <- data$best_models
  
  # List to store samples for each model
  lst_haz <- vector("list", length = nrow(dat))
  names(lst_haz) <- dat$ctgov
  
  for(i in seq_along(dat2$nested_params)) {
    
    # Get names, parameters and survival times
    name_loop <- dat$ctgov[i]
    params_tbl <- dat2$nested_params[[i]]
    surv_times <- dat$t2e[[i]]
    
    for(j in seq_along(dat2$nested_params[[i]]$params)) {
      
      # Get distribution and sequence times
      dist <- paste0(dat2$nested_params[[i]]$dist[[j]])
      time <- round(seq(min(surv_times$time), max(surv_times$time)))
      
      # Conditionally estimate hazard
      if(dist == "exp") {
        
        rate <- dat2$nested_params[[i]]$params[[j]]
        haz <- hexp(x = time, rate = exp(rate))
        
        lst_haz[[name_loop]][[dist]] <- data.frame(
          ctgov = dat$ctgov[i],
          dist = dist,
          haz_est = haz,
          time = time
        )
        
      } else if (dist == "gengamma") {
        
        mu <- as.data.frame(dat2$nested_params[[i]]$params[[j]])$mu
        sigma <- as.data.frame(dat2$nested_params[[i]]$params[[j]])$sigma
        q <- as.data.frame(dat2$nested_params[[i]]$params[[j]])$Q
        
        haz <- hgengamma(
          x = time, 
          mu = mu,
          sigma = exp(sigma),
          Q = q
        )
        
        lst_haz[[name_loop]][[dist]] <- data.frame(
          ctgov = dat$ctgov[i],
          dist = dist,
          haz_est = haz,
          time = time
        )
        
      } else if (dist == "gompertz") {
        
        shape <- as.data.frame(dat2$nested_params[[i]]$params[[j]])$shape
        rate <- as.data.frame(dat2$nested_params[[i]]$params[[j]])$rate
        
        haz <- hgompertz(
          x = time, 
          shape = shape,
          rate = exp(rate)
        )
        
        lst_haz[[name_loop]][[dist]] <- data.frame(
          ctgov = dat$ctgov[i],
          dist = dist,
          haz_est = haz,
          time = time
        )
        
      } else if (dist == "llogis") {
        
        shape <- as.data.frame(dat2$nested_params[[i]]$params[[j]])$shape
        scale <- as.data.frame(dat2$nested_params[[i]]$params[[j]])$scale
        
        haz <- hllogis(
          x = time, 
          shape = exp(shape),
          scale = exp(scale)
        )
        
        lst_haz[[name_loop]][[dist]] <- data.frame(
          ctgov = dat$ctgov[i],
          dist = dist,
          haz_est = haz,
          time = time
        )
        
      } else if (dist == "lnorm") {
        
        meanlog <- as.data.frame(dat2$nested_params[[i]]$params[[j]])$meanlog
        sdlog <- as.data.frame(dat2$nested_params[[i]]$params[[j]])$sdlog
        
        haz <- hlnorm(
          x = time, 
          meanlog = meanlog,
          sdlog = exp(sdlog)
        )
        
        lst_haz[[name_loop]][[dist]] <- data.frame(
          ctgov = dat$ctgov[i],
          dist = dist,
          haz_est = haz,
          time = time
        )
        
      } else if (dist == "weibull") {
        
        scale <- as.data.frame(dat2$nested_params[[i]]$params[[j]])$scale
        shape <- as.data.frame(dat2$nested_params[[i]]$params[[j]])$shape
        
        haz <- hweibull(
          x = time, 
          shape = exp(shape),
          scale = exp(scale)
        )
        
        lst_haz[[name_loop]][[dist]] <- data.frame(
          ctgov = dat$ctgov[i],
          dist = dist,
          haz_est = haz,
          time = time
        )
      }
    }
    
    lst_haz[[i]] <- bind_rows(lst_haz[[i]])
    
  }
  
  df_haz <- bind_rows(lst_haz) %>% arrange(ctgov)
  
  return(df_haz)
  
}








