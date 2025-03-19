
### Sample parameters and generate parametric estimates ------------------------

# Load packages, functions, etc.
source("Scripts/00a_packages_functions_coa.R")

# Import datasets
import_covariance <- read_csv("../Analysis/Vivli/Models_coa/covariance.csv")
import_fit <- read_csv("../Analysis/Vivli/Models_coa/fit.csv")
import_parameters <- read_csv("../Analysis/Vivli/Models_coa/parameters.csv")
import_t2e <- read_csv("../Analysis/Vivli/Models_coa/t2e.csv")

# Get transformed datasets
data <- get_data()

# Check Kaplan-Meier curves
plot_t2e <- data$df_t2e %>% 
  ggplot(
    aes(
      x = time,
      y = 1 - estimate,
      colour = cause
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
  "../Analysis/Output/kaplan_meier.png",
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
    paste0("../Analysis/Data/coa/", name, ".rds")
  )
}

### Plot parametric and KM estimates -------------------------------------------

# Load packages, functions, etc.
source("Scripts/00a_packages_functions_coa.R")

# Import datasets
data <- read_rds("../Analysis/Data/coa/data.rds")
sample_parameters <- read_rds("../Analysis/Data/coa/sample_parameters.rds")
sample_estimates <- read_rds("../Analysis/Data/coa/sample_estimates.rds")
plots_coord_cartesian <- read_csv("../Analysis/Data/coa/plots_for_coord_cartesian.csv")
import_flow <- read_csv(
  "../Analysis/Vivli/Summary_statistics/participant_flow_formatted.csv"
)

# Check plots where coord_cartesian needs applied
plot_t2e <- data$df_t2e %>% 
  mutate(ctgov_cause = paste0(ctgov, "_", cause)) %>% 
  filter(ctgov_cause %in% plots_coord_cartesian$ctgov_cause) %>% 
  ggplot(
    aes(
      x = time,
      y = 1 - estimate,
      colour = cause
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
  "../Analysis/Output/kaplan_meier_coord_targets.png",
  plot_t2e,
  dpi = 300,
  width = 18,
  height = 12
)

# Get plots of Kaplan-Meier and parametric estimates
plot_estimates <- get_plots()

# Gt table of plots
gt_plots <- get_gt(type = "all")

# Save gt of all trials
gtsave(gt_plots, "../Analysis/Output/gt_all.html")

# Split trials based on number of causes
# Part 1: Trials where at least 5/7 causes occurred (N=19)
# Part 2: Remaining trials (N=72)

# Select trials
split_ref <- unnest(data$nested_t2e, cols = c("t2e")) %>% 
  inner_join(import_flow %>% dplyr::select(ctgov)) %>% 
  group_by(ctgov, cause) %>% 
  reframe() %>% 
  group_by(ctgov) %>% 
  reframe(n_causes = n())

split_main <- split_ref %>% filter(n_causes >= 5)
split_supp <- split_ref %>% filter(n_causes < 5)

# Get main trial plots
# Order list names
plots_main <- list()

for (ctgov in split_main$ctgov) {
  match <- grep(ctgov, names(plot_estimates), value = TRUE)
  for(name in match) {
    plots_main[[name]] <- plot_estimates[[name]]
  }
}

var_names_main <- data.frame(nms = names(plots_main)) %>% 
  separate(
    nms,
    into = c("ctgov", "cause", "dist"),
    sep = "_"
  ) %>% 
  mutate(
    ctgov_cause = paste0(ctgov, "_", cause, "_", dist)
  ) %>% 
  distinct()

var_names_main_order <- import_flow %>% 
  dplyr::select(ctgov, condition) %>%
  inner_join(var_names_main) %>% 
  arrange(condition, ctgov)

var_main <- unique(var_names_main_order$ctgov_cause)

plots_main_ordered <- plots_main[match(var_main, names(plots_main))]

# Get supp trial plots
# Order list names
plots_supp <- list()

for (ctgov in split_supp$ctgov) {
  match <- grep(ctgov, names(plot_estimates), value = TRUE)
  for(name in match) {
    plots_supp[[name]] <- plot_estimates[[name]]
  }
}

var_names_supp <- data.frame(nms = names(plots_supp)) %>% 
  separate(
    nms,
    into = c("ctgov", "cause", "dist"),
    sep = "_"
  ) %>% 
  mutate(
    ctgov_cause = paste0(ctgov, "_", cause, "_", dist)
  ) %>% 
  distinct()

var_names_supp_order <- import_flow %>% 
  dplyr::select(ctgov, condition) %>%
  inner_join(var_names_supp) %>% 
  arrange(condition, ctgov)

var_supp <- unique(var_names_supp_order$ctgov_cause)

plots_supp_ordered <- plots_supp[match(var_supp, names(plots_supp))]

# Get gt tables of split
gt_main <- get_gt(type = "main")
gt_supp <- get_gt(type = "supp")

# Save splits
gtsave(gt_main, "../Analysis/Output/gt/gt_main.html")
gtsave(gt_supp, "../Analysis/Output/gt/gt_supp.html")

# Split further for supplementary figure groups
supp_p1 <- plots_supp_ordered[1:98] #20 - #33
supp_p2 <- plots_supp_ordered[99:196] #34 - #47
supp_p3 <- plots_supp_ordered[197:294] #48 - #62
supp_p4 <- plots_supp_ordered[295:392] #63 - #75
supp_p5 <- plots_supp_ordered[393:497] #76 - #90

supp_p1 <- get_gt(type = "supp_parts", part = 1)
gtsave(supp_p1, "../Analysis/Output/gt/gt_supp_pt1.html")

supp_p2 <- get_gt(type = "supp_parts", part = 2)
gtsave(supp_p2, "../Analysis/Output/gt/gt_supp_pt2.html")

supp_p3 <- get_gt(type = "supp_parts", part = 3)
gtsave(supp_p3, "../Analysis/Output/gt/gt_supp_pt3.html")

supp_p4 <- get_gt(type = "supp_parts", part = 4)
gtsave(supp_p4, "../Analysis/Output/gt/gt_supp_pt4.html")

supp_p5 <- get_gt(type = "supp_parts", part = 5)
gtsave(supp_p5, "../Analysis/Output/gt/gt_supp_pt5.html")

# Plot hazard rates ------------------------------------------------------------

# Load packages, functions, etc.
source("Scripts/00a_packages_functions_coa.R")

df_hazard <- get_haz()
df_hazard <- df_hazard %>% 
  filter(
    dist != "exp",
    !(ctgov %in% c("NCT01131676", "NCT00274573"))
  )

# Plot hazard rate
plot_haz <- df_hazard %>% 
  group_by(ctgov, cause) %>% 
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
  facet_wrap(~cause, scales = "free", nrow = 2) +
  theme_bw() +
  theme(
    text = element_text(size = 6)
  )

ggsave(
  "../Analysis/Output/plots/jce_hazard_causes.png",
  plot_haz,
  dpi = 300,
  width = 2244/300,
  height = 1683/300,
  units = "in"
)

# Plot hazard rate for tally count
plot_haz2 <- df_hazard %>% 
  group_by(ctgov, cause) %>% 
  mutate(
    est = haz_est/max(haz_est)
  ) %>% 
  ggplot(
    aes(
      x = time,
      y = est,
      colour = cause
    )
  ) +
  geom_line(linewidth = 1) +
  labs(
    x = "Trial duration in days",
    y = "Hazard estimates as a proportion of highest point estimate"
  ) +
  facet_wrap(~ctgov, scales = "free") +
  theme_bw()

ggsave(
  "../Analysis/Output/plots/hazard_tally.png",
  plot_haz2,
  dpi = 300,
  width = 18,
  height = 12
)



