
### Sample parameters and generate parametric estimates ------------------------

# Load packages, functions, etc.
source("Scripts/00b_packages_functions_total.R")

# Import datasets
import_covariance <- read_csv("../Analysis/Vivli/Models_ta/covariance.csv")
import_fit <- read_csv("../Analysis/Vivli/Models_ta/fit.csv")
import_parameters <- read_csv("../Analysis/Vivli/Models_ta/parameters.csv")
import_t2e <- read_csv("../Analysis/Vivli/Models_ta/t2e.csv")

# Get transformed datasets
data <- get_data()

# Check Kaplan-Meier curves
plot_t2e <- data$df_t2e %>% 
  ggplot(
    aes(
      x = time,
      y = 1 - estimate
    )
  ) +
  geom_step(linewidth = 1) +
  facet_wrap(
    ~ctgov, 
    scales = "free"
  ) +
  labs(
    x = "Time in trial (Days)",
    y = "Probability of attrition occurring"
  )

# Save plot
ggsave(
  "../Analysis/Output/plots/kaplan_meier_ta.png",
  plot_t2e,
  dpi = 300,
  width = 18,
  height = 12
)

# Sample parameters from a multivariate normal distribution 
sample_parameters <- get_samples()

# Generate parametric estimates using sampled parameters
sample_estimates <- get_estimates()

# Save outputs
save_names <- c(
  "data", "sample_parameters", "sample_estimates"
)
save_lst <- lapply(save_names, function(i) get(i))
names(save_lst) <- save_names

for (i in seq_along(save_lst)) {
  name <- names(save_lst[i])
  msg <- paste0("Saving ", name)
  print(msg)
  write_rds(
    save_lst[[i]],
    paste0("../Analysis/Data/ta/", name, ".rds")
  )
}

### Plot parametric and KM estimates -------------------------------------------

# Load packages, functions, etc.
source("Scripts/00b_packages_functions_total.R")

# Import datasets
data <- read_rds("../Analysis/Data/ta/data.rds")
sample_parameters <- read_rds("../Analysis/Data/ta/sample_parameters.rds")
sample_estimates <- read_rds("../Analysis/Data/ta/sample_estimates.rds")
import_flow <- read_csv("../Analysis/Vivli/Summary_statistics/participant_flow_formatted.csv")

# Bind estimates for plotting
sample_estimates_bind <- lapply(seq_along(sample_estimates), function(i) {
  name <- names(sample_estimates[i])
  lapply(seq_along(sample_estimates[[i]]), function(j) {
    name2 <- names(sample_estimates[[i]][j])
    bind_rows(
      sample_estimates[[i]][[j]]
    ) %>% 
      mutate(
        ctgov = name,
        dist = name2
      )
  })
})

sample_estimates_plot <- bind_rows(sample_estimates_bind) %>% 
  filter(
    !is.na(estimate),
    !(ctgov %in% c("NCT01131676", "NCT00274573"))
  )

# Get T2E plots
total_plot <- sample_estimates_plot %>% 
  # mutate(
  #   new_id = dense_rank(ctgov),
  #   new_id = if_else(new_id < 10L, paste0("0", new_id), as.character(new_id)),
  #   new_id = paste0("Trial ", new_id),
  #   new_id = factor(new_id)
  # ) |> 
  ggplot(
    aes(
      x = time,
      y = round((1 - estimate) * 100, 1),
      group = iter
    )
  ) +
  geom_step(
    linewidth = 0.5,
    colour = "black"
  ) +
  geom_line(
    aes(
      x = time, 
      y = round((1 - fit_est) * 100, 1)
    ), 
    linewidth = 0.5, 
    alpha = 0.02,
    colour = case_when(
      sample_estimates_plot$dist == "exp" ~ "skyblue",
      #sample_estimates_plot$dist == "gengamma" ~ "limegreen",
      sample_estimates_plot$dist == "gompertz" ~"orange",
      sample_estimates_plot$dist == "lnorm" ~ "purple",
      sample_estimates_plot$dist == "llogis" ~ "deeppink",
      sample_estimates_plot$dist == "weibull" ~ "red",
      TRUE ~ "black"
    )
  ) +
  scale_x_continuous(n.breaks = 4) +
  scale_y_continuous(n.breaks = 4) +
  labs(
    x = "Time in days",
    y = "Cumulative incidence of attrition (%)"
  ) +
  facet_wrap(
    ~ ctgov,
    scales = "free"
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 6)
  )

ggsave(
  "../Analysis/Output/plots/jce_t2e_ta.png",
  total_plot,
  dpi = 300,
  width = 2244/300,
  height = 1683/300,
  units = "in"
)  

single_plot <- sample_estimates_plot %>% 
  mutate(
    new_id = dense_rank(ctgov),
    new_id = if_else(new_id < 10L, paste0("0", new_id), as.character(new_id)),
    new_id = paste0("Trial ", new_id),
    new_id = factor(new_id)
  ) |> 
  filter(new_id == "Trial 26") |> 
  ggplot(
    aes(
      x = time,
      y = round((1 - estimate) * 100, 1),
      group = iter
    )
  ) +
  geom_step(
    linewidth = 0.5,
    colour = "black"
  ) +
  geom_line(
    aes(
      x = time, 
      y = round((1 - fit_est) * 100, 1)
    ), 
    linewidth = 0.5, 
    alpha = 0.05,
    colour = "orange"
  ) +
  scale_x_continuous(n.breaks = 6) +
  scale_y_continuous(n.breaks = 5) +
  labs(
    x = "Time in days",
    y = "Cumulative incidence of attrition (%)"
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 9)
  )

ggsave(
  "../Analysis/Output/ictmc/example.png",
  single_plot,
  dpi = 300,
  width = 4,
  height = 4
)

ggsave(
  "../Analysis/Output/plots/t2e_ta.png",
  total_plot,
  dpi = 300,
  width = 14,
  height = 8
)

# Plot hazard rates ------------------------------------------------------------

# Load packages, functions, etc.
source("Scripts/00b_packages_functions_total.R")

df_hazard <- get_haz()
df_hazard <- df_hazard %>% 
  filter(
    dist != "exp", 
    !(ctgov %in% c("NCT01131676", "NCT00274573"))
  )

# Plot hazard rate
plot_haz <- df_hazard %>% 
  # mutate(
  #   new_id = dense_rank(ctgov),
  #   new_id = if_else(new_id < 10L, paste0("0", new_id), as.character(new_id)),
  #   new_id = paste0("Trial ", new_id),
  #   new_id = factor(new_id)
  # ) |> 
  group_by(ctgov) %>% 
  mutate(
    est = haz_est/max(haz_est)
  ) %>% 
  ggplot(
    aes(
      x = time,
      y = est,
      group = ctgov,
      colour = dist
    )
  ) +
  geom_line(
    linewidth = 0.5,
    colour = case_when(
      #df_hazard$dist == "exp" ~ "skyblue",
      #df_hazard$dist == "gengamma" ~ "limegreen",
      df_hazard$dist == "gompertz" ~"orange",
      df_hazard$dist == "lnorm" ~ "purple",
      df_hazard$dist == "llogis" ~ "deeppink",
      df_hazard$dist == "weibull" ~ "red",
      TRUE ~ "black"
    )
  ) +
  scale_x_continuous(n.breaks = 4) +
  scale_y_continuous(n.breaks = 4) +
  labs(
    x = "Time in days",
    y = "Scaled hazard rate"
  ) +
  facet_wrap(~ctgov, scales = "free") +
  theme_bw() +
  theme(
    text = element_text(size = 6)
  )

ggsave(
  "../Analysis/Output/plots/jce_hazard_ta.png",
  plot_haz,
  dpi = 300,
  width = 2244/300,
  height = 1683/300,
  units = "in"
)

single_plot <- df_hazard %>% 
  mutate(
    new_id = dense_rank(ctgov),
    new_id = if_else(new_id < 10L, paste0("0", new_id), as.character(new_id)),
    new_id = paste0("Trial ", new_id),
    new_id = factor(new_id)
  ) |> 
  group_by(new_id) %>% 
  mutate(
    est = haz_est/max(haz_est)
  ) %>% 
  ungroup() |> 
  filter(new_id == "Trial 22") |> 
  ggplot(
    aes(
      x = time,
      y = est,
    )
  ) +
  geom_line(
    linewidth = 1,
    colour = "orange"
  ) +
  labs(
    x = "Time in days",
    y = "Scaled hazard rate"
  ) +
  theme_bw()

ggsave(
  "../Analysis/Output/ictmc/example_haz.png",
  single_plot,
  dpi = 300,
  width = 4,
  height = 4
)
