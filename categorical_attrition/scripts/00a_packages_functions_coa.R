
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
    group_by(ctgov, cause) %>% 
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
      ),
      ctgov_cause = paste0(ctgov, "_", cause)
    )
  
  # Nest Kaplan-Meier data
  t2e <- t2e_transform %>% 
    group_by(ctgov, cause) %>% 
    nest(.key = "t2e") %>% 
    ungroup()
  
  # Join parameters and fit metrics
  params_fit <- dat2 %>% left_join(dat3) %>% arrange(ctgov, cause)
  
  print("Getting best models")
  
  # Select the best 3 models per cause and trial based on AIC
  # Many cases where 4 are selected because the AIC values were the same between
  # two or more
  get_best_models <- params_fit %>% 
    dplyr::select(ctgov, cause, dist, aic) %>% 
    mutate(ctgov_cause = paste0(ctgov, "_", cause)) %>% 
    filter(dist != "gengamma") %>% 
    dplyr::select(-ctgov_cause) %>% 
    distinct() %>% 
    group_by(ctgov, cause) %>% 
    slice_min(aic, n = 1) %>% 
    slice_head(n = 1) %>% 
    ungroup()
  
  print("Nesting best model parameters")
  
  # Nest parameters for best performing models
  best_models <- params_fit %>% 
    inner_join(
      get_best_models %>% 
        dplyr::select(ctgov, cause, dist)
    ) %>% 
    dplyr::select(-aic, -loglik) %>% 
    nest(
      .key = "nested_params",
      .by = c("ctgov", "cause")
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
    arrange(ctgov, cause)
  
  print("Nesting variance-covariance")
  
  # Nest variance-covariance matrices
  vcov <- import_covariance %>% 
    inner_join(
      get_best_models %>% 
        dplyr::select(ctgov, cause, dist)
    ) %>% 
    nest(
      .key = "distr_cov",
      .by = c("ctgov", "cause")
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
    arrange(ctgov, cause)
  
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
  names(lst_samples) <- paste0(dat$ctgov, "_", dat$cause)
  
  for(i in seq_along(dat2$nested_params)) {
    
    name_loop <- paste0(dat$ctgov[i], "_", dat$cause[i])
    
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
  names(lst_estimates) <- paste0(dat$ctgov, "_", dat$cause)
  
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

# Get plots of Kaplan-Meier and parametric estimates ---------------------------

get_plots <- function() {
  
  dat <- sample_estimates[!grepl("NCT01131676", names(sample_estimates))]
  dat2 <- data$nested_t2e %>% filter(ctgov != "NCT01131676")
  dat3 <- plots_coord_cartesian %>% filter(!grepl("NCT01131676", ctgov_cause))
  
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

# Gt table of plots ------------------------------------------------------------

get_gt <- function(type = c("all", "main", "supp", "supp_parts"), part) {
  
  # Conditionally assign data
  dat <- data$nested_t2e %>% 
    filter(ctgov != "NCT01131676")
  
  if (type == "all") {
    dat2 <- plot_estimates
  } else if (type == "main") {
    dat2 <- plots_main_ordered
  } else if (type == "supp") {
    dat2 <- plots_supp_ordered
  } else if (type == "supp_parts" & part == 1) {
    dat2 <- supp_p1
  } else if (type == "supp_parts" & part == 2) {
    dat2 <- supp_p2
  } else if (type == "supp_parts" & part == 3) {
    dat2 <- supp_p3
  } else if (type == "supp_parts" & part == 4) {
    dat2 <- supp_p4
  } else if (type == "supp_parts" & part == 5) {
    dat2 <- supp_p5
  }
  
  # Get ctgov/conditions ordered properly
  names_plots <- data.frame(
    names = names(dat2)
  ) %>% 
    separate(
      names, 
      into = c("ctgov", "cause", "dist"), 
      sep = "_"
    )
  
  import_cond <- read_csv(
    "../Analysis/Vivli/Summary_statistics/participant_flow_formatted.csv"
  )
  
  ctgov <- dat %>% 
    filter(ctgov %in% names_plots$ctgov) %>% 
    inner_join(
      import_cond %>% 
        dplyr::select(ctgov, condition)
    ) %>% 
    dplyr::select(ctgov, condition) %>% 
    distinct()
  
  order_names_conditions <- ctgov %>% arrange(condition)
  
  gt <- dplyr::tibble(
    id = order_names_conditions$ctgov,
    cond = order_names_conditions$condition,
    ae_surv = NA,
    loe_surv = NA,
    l2f_surv = NA,
    o_surv = NA,
    pi_surv = NA,
    pd_surv = NA,
    vw_surv = NA
  ) %>%
    gt() %>%
    tab_style(  # Format value alignment and appearance for column labels
      style = cell_text(
        weight = "bold",
        align = "center",
        v_align = "middle"
      ),
      locations = cells_column_labels()
    ) %>% 
    cols_label(  # Add column labels for format of each cell
      id = "NCT ID",
      cond = "Condition studied",
      ae_surv = "Adverse Event",
      loe_surv = "Lack of Efficacy",
      l2f_surv = "Lost to Follow-up",
      pi_surv = "PI/Sponsor Decision",
      pd_surv = "Protocol Violation",
      o_surv = "Other",
      vw_surv = "Voluntary Withdrawal"
    ) %>%
    tab_options(  # Reduce row height
      data_row.padding = px(1)
    ) %>%
    tab_style(  # Add row dividers to title row
      style = cell_borders(
        sides = c("top", "bottom", "left", "right"),
        color = "black",
        weight = px(2)
      ),
      locations = cells_column_labels()
    ) %>% 
    tab_style(  # Add column dividers to title row
      style = cell_borders(
        sides = c("top", "bottom", "left", "right"),
        color = "black",
        weight = px(2)
      ),
      locations = cells_column_labels()
    ) %>%
    tab_style(  # Add column dividers to rest of rows
      style = cell_borders(
        sides = c("top", "bottom", "left", "right"),
        color = "black",
        weight = px(2)
      ),
      locations = cells_body(
        columns = c(
          id,
          cond,
          ae_surv,
          loe_surv,
          l2f_surv,
          o_surv,
          pi_surv,
          pd_surv,
          vw_surv
        )
      )
    ) %>% 
    text_transform(  # Adverse event
      locations = cells_body(columns = ae_surv),
      fn = function(x) {
        dat2[grep("Adverse", names(dat2))] %>%
          ggplot_image(height = px(30), aspect_ratio = 2)
      }
    ) %>%
    text_transform(  # Lack of Efficacy
      locations = cells_body(columns = loe_surv),
      fn = function(x) {
        dat2[grep("Lack", names(dat2))] %>%
          ggplot_image(height = px(30), aspect_ratio = 2)
      }
    ) %>%
    text_transform(  # Lost to Follow-up
      locations = cells_body(columns = l2f_surv),
      fn = function(x) {
        dat2[grep("Lost", names(dat2))] %>%
          ggplot_image(height = px(30), aspect_ratio = 2)
      }
    ) %>%
    text_transform(  # Other
      locations = cells_body(columns = o_surv),
      fn = function(x) {
        dat2[grep("Other", names(dat2))] %>%
          ggplot_image(height = px(30), aspect_ratio = 2)
      }
    ) %>%
    text_transform(  # PI/Sponsor Decision
      locations = cells_body(columns = pi_surv),
      fn = function(x) {
        dat2[grep("PI", names(dat2))] %>%
          ggplot_image(height = px(30), aspect_ratio = 2)
      }
    ) %>%
    text_transform(  # Protocol Deviation
      locations = cells_body(columns = pd_surv),
      fn = function(x) {
        dat2[grep("Protocol", names(dat2))] %>%
          ggplot_image(height = px(30), aspect_ratio = 2)
      }
    ) %>%
    text_transform(  # Voluntary Withdrawal
      locations = cells_body(columns = vw_surv),
      fn = function(x) {
        dat2[grep("Voluntary", names(dat2))] %>%
          ggplot_image(height = px(30), aspect_ratio = 2)
      }
    )
  
  return(gt)
  
}

# Get hazard estimates from model parameters -----------------------------------

get_haz <- function() {
  
  # Conditionally import data
  dat <- data$nested_t2e
  dat2 <- data$best_models
  
  # List to store samples for each model
  lst_haz <- vector("list", length = nrow(dat))
  names(lst_haz) <- paste0(dat$ctgov, "_", dat$cause)
  
  for(i in seq_along(dat2$nested_params)) {
    
    # Get names, parameters and survival times
    name_loop <- paste0(dat$ctgov[i], "_", dat$cause[i])
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
          cause = dat$cause[i],
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
          cause = dat$cause[i],
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
          cause = dat$cause[i],
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
          cause = dat$cause[i],
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
          cause = dat$cause[i],
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








