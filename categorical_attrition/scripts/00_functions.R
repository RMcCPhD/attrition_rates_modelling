
# Sample parameters from a multivariate normal distribution --------------------

get_samples <- function() {
  
  # Define data
  dat <- lst_prep$cind
  dat2 <- lst_prep$params
  dat3 <- lst_prep$vcov
  
  # List to store samples for each model
  lst_samples <- vector("list", length = nrow(dat))
  names(lst_samples) <- paste0(dat$ctgov, "_", dat$cause)
  
  for(i in seq_along(dat2$params)) {
    
    name_loop <- paste0(dat$ctgov[i], "_", dat$cause[i])
    
    for(j in seq_along(dat2$params[[i]]$params)) {
      
      dist <- paste0(dat2$params[[i]]$dist[[j]])
      surv <- dat$cind[[i]]
      
      n = 100 # Number of samples
      mu = dat2$params[[i]]$params[[j]]  # Parameter estimates
      sigma = dat3$vcov[[i]]$vcov[[j]]  # Variance-covariance
      
      # Sample parameters with mvrnorm and catch any errors
      tryCatch({
          lst_samples[[i]][[dist]] <- mvrnorm(n = n, mu = mu, Sigma = sigma)
        },
        error = function(e) {
          message(paste0("Error sampling ", name_loop, " ", dist, " Error: ", e))
          return(NULL)
        })
    }
  }
  
  return(lst_samples)
  
}

# Generate cumulative incidence estimates using sampled parameters -------------

get_cind <- function() {
  dat <- lst_prep$cind
  dat2 <- sample_parameters
  
  lst_estimates <- vector("list", length = nrow(dat))
  names(lst_estimates) <- paste0(dat$ctgov, "_", dat$cause)
  
  for(i in seq_along(dat2)) {
    
    cind <- dat$cind[[i]]
    name_loop <- names(dat2[i])
    print(name_loop)
    
    time <- round(seq(min(cind$time), max(cind$time)))
    
    cind_by_time <- cind
    
    for(j in seq_along(dat2[[i]])) {
      dist <- names(dat2[[i]][j])
      print(dist)
      
      params_df <- as.data.frame(dat2[[i]][[j]])
      dist_results <- vector("list", nrow(params_df))
      
      for(k in seq_len(nrow(params_df))) {
        name_k <- paste0("iter_", k)
        
        fit_est <- switch(
          dist,
          exp = 1 - pexp(time, rate = exp(params_df[k, 1])),
          gompertz = 1 - pgompertz(time, rate = exp(params_df[k, "rate"]), shape = params_df[k, "shape"]),
          llogis = 1 - pllogis(time, scale = exp(params_df[k, "scale"]), shape = exp(params_df[k, "shape"])),
          lnorm = 1 - plnorm(time, meanlog = params_df[k, "meanlog"], sdlog = exp(params_df[k, "sdlog"])),
          weibull = 1 - pweibull(time, scale = exp(params_df[k, "scale"]), shape = exp(params_df[k, "shape"])),
          stop(paste("Unknown distribution:", dist))
        )
        
        dist_results[[name_k]] <- data.frame(
          iter = k,
          fit_est = fit_est,
          time = time
        )
      }
      
      # Bind results for this distribution and join once
      all_results <- dplyr::bind_rows(dist_results)
      all_results <- dplyr::left_join(all_results, cind_by_time, by = "time")
      lst_estimates[[name_loop]][[dist]] <- split(all_results, all_results$iter)
    }
  }
  
  return(lst_estimates)
}

# Get hazard estimates from model parameters -----------------------------------

get_hazards <- function() {
  
  # Conditionally import data
  dat <- lst_prep$cind
  dat2 <- lst_prep$params
  
  # List to store samples for each model
  lst_haz <- vector("list", length = nrow(dat))
  names(lst_haz) <- paste0(dat$ctgov, "_", dat$cause)
  
  for(i in seq_along(dat2$params)) {
    
    # Get names, parameters and survival times
    name_loop <- paste0(dat$ctgov[i], "_", dat$cause[i])
    params_tbl <- dat2$params[[i]]
    surv_times <- dat$cind[[i]]
    
    for(j in seq_along(dat2$params[[i]]$params)) {
      
      # Get distribution and sequence times
      dist <- paste0(dat2$params[[i]]$dist[[j]])
      time <- round(seq(min(surv_times$time), max(surv_times$time)))
      
      # Conditionally estimate hazard
      if(dist == "exp") {
        
        rate <- dat2$params[[i]]$params[[j]]
        haz <- hexp(x = time, rate = exp(rate))
        
        lst_haz[[name_loop]][[dist]] <- data.frame(
          ctgov = dat$ctgov[i],
          cause = dat$cause[i],
          dist = dist,
          haz_est = haz,
          time = time
        )
        
      } else if (dist == "gompertz") {
        
        shape <- as.data.frame(dat2$params[[i]]$params[[j]])$shape
        rate <- as.data.frame(dat2$params[[i]]$params[[j]])$rate
        
        haz <- hgompertz(
          x = time, 
          shape = shape,
          rate = exp(rate)
        )
        
        lst_haz[[name_loop]][[dist]] <- data.frame(
          ctgov = dat$ctgov[i],
          cause = dat$cause[i],
          dist = dist,
          haz_est = haz,
          time = time
        )
        
      } else if (dist == "llogis") {
        
        shape <- as.data.frame(dat2$params[[i]]$params[[j]])$shape
        scale <- as.data.frame(dat2$params[[i]]$params[[j]])$scale
        
        haz <- hllogis(
          x = time, 
          shape = exp(shape),
          scale = exp(scale)
        )
        
        lst_haz[[name_loop]][[dist]] <- data.frame(
          ctgov = dat$ctgov[i],
          cause = dat$cause[i],
          dist = dist,
          haz_est = haz,
          time = time
        )
        
      } else if (dist == "lnorm") {
        
        meanlog <- as.data.frame(dat2$params[[i]]$params[[j]])$meanlog
        sdlog <- as.data.frame(dat2$params[[i]]$params[[j]])$sdlog
        
        haz <- hlnorm(
          x = time, 
          meanlog = meanlog,
          sdlog = exp(sdlog)
        )
        
        lst_haz[[name_loop]][[dist]] <- data.frame(
          ctgov = dat$ctgov[i],
          cause = dat$cause[i],
          dist = dist,
          haz_est = haz,
          time = time
        )
        
      } else if (dist == "weibull") {
        
        scale <- as.data.frame(dat2$params[[i]]$params[[j]])$scale
        shape <- as.data.frame(dat2$params[[i]]$params[[j]])$shape
        
        haz <- hweibull(
          x = time, 
          shape = exp(shape),
          scale = exp(scale)
        )
        
        lst_haz[[name_loop]][[dist]] <- data.frame(
          ctgov = dat$ctgov[i],
          cause = dat$cause[i],
          dist = dist,
          haz_est = haz,
          time = time
        )
      }
    }
    
    lst_haz[[i]] <- bind_rows(lst_haz[[i]])
    
  }
  
  df_haz <- bind_rows(lst_haz) %>% arrange(ctgov, cause)
  
  return(df_haz)
  
}








