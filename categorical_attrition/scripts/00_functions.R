
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

# Get cumulative incidence plots for arranging with ggtable --------------------

get_plots <- function() {
  
  dat <- sample_estimates
  dat2 <- lst_prep$cind
  dat3 <- plots_coord_cartesian
  
  # List to store plots for each dist per trial/cause
  lst_plots <- vector("list")
  
  for(i in seq_along(dat)) {
    
    name_loop <- names(dat[i])
    print(name_loop)
    
    for(j in seq_along(dat[[i]])) {
      
      dist <- names(dat[[i]][j])
      plot_dat <- bind_rows(dat[[i]][[j]])
      lst_name <- paste0(name_loop, "_", dist)
      plot_name <- gsub("^(.*?)_(.*?)_.*$", "\\1 - \\2", lst_name)
      
      print(lst_name)
      
      lst_plots[[lst_name]] <- plot_dat %>% 
        filter(
          !is.na(estimate)
        ) %>% 
        ggplot(
          aes(
            x = time, 
            y = 1 - estimate, 
            group = iter
          )
        ) +
        geom_step(
          linewidth = 12, 
          colour = "black"
        ) +
        geom_line(
          aes(
            x = time, 
            y = 1 - fit_est
          ), 
          linewidth = 12, 
          colour = case_when(
            grepl("exp", lst_name) ~ "skyblue",
            grepl("gengamma", lst_name) ~ "limegreen",
            grepl("gompertz", lst_name) ~"orange",
            grepl("lnorm", lst_name) ~ "purple",
            grepl("llogis", lst_name) ~ "salmon",
            grepl("weibull", lst_name) ~ "red",
            TRUE ~ "black"
          ), 
          alpha = 0.025
        ) +
        scale_x_continuous(
          limits = c(min(plot_dat$time), max(plot_dat$time)),
          n.breaks = 6
        ) +
        guides(alpha = "none") +
        theme_bw() +
        theme_void()
    }
  }
  
  # Get list names for each applications of coord_cartesian
  names_01 <- dat3 %>% filter(apply_01 == 1) %>% dplyr::select(ctgov_cause)
  names_02 <- dat3 %>% filter(apply_02 == 1) %>% dplyr::select(ctgov_cause)
  names_03 <- dat3 %>% filter(apply_03 == 1) %>% dplyr::select(ctgov_cause)
  names_05 <- dat3 %>% filter(apply_05 == 1) %>% dplyr::select(ctgov_cause)
  names_06 <- dat3 %>% filter(apply_06 == 1) %>% dplyr::select(ctgov_cause)
  names_10 <- dat3 %>% filter(apply_10 == 1) %>% dplyr::select(ctgov_cause)
  names_15 <- dat3 %>% filter(apply_15 == 1) %>% dplyr::select(ctgov_cause)
  
  # Apply coord_cartesian where estimates diverge considerably
  lst_plots_coord <- lapply(names(lst_plots), function(plot_name) {
    if (any(grepl(paste(names_01$ctgov_cause, collapse = "|"), plot_name))) {
      lst_plots[[plot_name]] + coord_cartesian(ylim = c(0, 0.01))
      
    } else if (any(grepl(paste(names_02$ctgov_cause, collapse = "|"), plot_name))) {
      lst_plots[[plot_name]] + coord_cartesian(ylim = c(0, 0.02))
      
    } else if (any(grepl(paste(names_03$ctgov_cause, collapse = "|"), plot_name))) {
      lst_plots[[plot_name]] + coord_cartesian(ylim = c(0, 0.03))
      
    } else if (any(grepl(paste(names_05$ctgov_cause, collapse = "|"), plot_name))) {
      lst_plots[[plot_name]] + coord_cartesian(ylim = c(0, 0.05))
      
    } else if (any(grepl(paste(names_06$ctgov_cause, collapse = "|"), plot_name))) {
      lst_plots[[plot_name]] + coord_cartesian(ylim = c(0, 0.06))
      
    } else if (any(grepl(paste(names_10$ctgov_cause, collapse = "|"), plot_name))) {
      lst_plots[[plot_name]] + coord_cartesian(ylim = c(0, 0.10))
      
    } else if (any(grepl(paste(names_15$ctgov_cause, collapse = "|"), plot_name))) {
      lst_plots[[plot_name]] + coord_cartesian(ylim = c(0, 0.15))
      
    } else {
      lst_plots[[plot_name]]
    }
  })
  
  names(lst_plots_coord) <- names(lst_plots)
  lst_plots <- lst_plots_coord
  
  # Create blank plots for where a cause did not occur
  ctgov <- unique(dat2$ctgov)
  cause <- unique(dat2$cause)
  
  combine <- expand.grid(
    ctgov = ctgov, 
    cause = cause
  ) %>%
    mutate(
      nms = paste0(ctgov, "_", cause)
    )
  
  blank_plots <- lapply(seq_along(combine$nms), function(i) {
    ggplot() +
      theme_void()
  })
  names(blank_plots) <- paste0(combine$nms)
  
  # Filter for blank plots that are not present in list of plots
  names_main <- gsub("^(.*?)_(.*?)_.*$", "\\1_\\2", names(lst_plots))
  blanks <- combine %>% filter(!nms %in% names_main)
  blank_plots_filter <- blank_plots[names(blank_plots) %in% blanks$nms]
  names(blank_plots_filter) <- paste0(blanks$nms, "_blank")
  
  # Join blanks to list of plots
  lst_plots <- c(lst_plots, blank_plots_filter)
  order_strings <- gsub("^(.*?)_(.*?)_.*$", "\\1_\\2", names(lst_plots))
  lst_plots <- lst_plots[order(order_strings)]
  
  # Order plots by condition, then ctgov
  lst_names <- data.frame(names = names(lst_plots)) %>% 
    separate(names, into = c("ctgov", "cause", "dist"), sep = "_")
  
  import_cond <- read_csv("../Analysis/Vivli/Summary_statistics/participant_flow_formatted.csv")
  
  ctgov <- dat2 %>% 
    inner_join(
      import_cond %>% 
        dplyr::select(ctgov, condition)
    ) %>% 
    dplyr::select(ctgov, condition) %>% 
    distinct()
  
  order_names_conditions <- ctgov %>% arrange(condition)
  
  names_conditions <- lst_names %>% 
    left_join(order_names_conditions) %>% 
    arrange(condition) %>% 
    mutate(names = paste0(ctgov, "_", cause, "_", dist))
  
  sort_names_conditions <- names_conditions$names
  
  lst_plots <- lst_plots[match(sort_names_conditions, names(lst_plots))]
  
  return(lst_plots)
  
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








