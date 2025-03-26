
source("config.R")
source("categorical_attrition/scripts/00_functions.R")
source("categorical_attrition/scripts/00_materials.R")

library(tidyverse)
library(flexsurv)
library(gt)
library(webshot2)

# Import data
imp_params <- read_csv("categorical_attrition/data/parameters.csv")
lst_prep <- read_rds("categorical_attrition/processed_data/prepared_data.rds")
imp_cind <- read_rds("categorical_attrition/processed_data/cind_sampled_df.rds")
imp_sum <- read_csv("vivli/data/summaries/participant_flow_formatted.csv")

# Create abbreviations for conditions
cond_short <- imp_sum %>% 
  group_by(condition) %>% 
  reframe() %>% 
  mutate(
    condition_short = case_match(
      condition,
      "Asthma" ~ "Asth",
      "Benign prostatic hyperplasia" ~ "BPH",
      "Hypertension" ~ "Hyp",
      "Osteoarthritis" ~ "Oarth",
      "Osteoporosis" ~ "Oporo",
      "Parkinson's disease" ~ "Park",
      "Pulmonary fibrosis" ~ "Pfib",
      "Restless leg syndrome" ~ "RLS",
      "Type 2 diabetes" ~ "T2DM",
      .default = condition
    )
  )

# Create index for grouping by condition and sample size
ind_cond_n <- imp_sum %>% 
  select(ctgov, participants, condition) %>% 
  rename(n = participants) %>% 
  left_join(cond_short) %>% 
  arrange(condition, desc(n)) %>% 
  group_by(condition) %>% 
  mutate(
    trial_num = row_number(),
    trial_num = if_else(trial_num < 10, paste0("0", trial_num), as.character(trial_num)),
    plot_id = paste0(condition_short, " ", trial_num, ", n=", n)
  ) %>% 
  ungroup() %>% 
  select(-trial_num)

# Join new id to cumulative incidence and hazard rate datasets
cind_plot_df <- ind_cond_n %>% left_join(imp_cind)

# Prepare dataset
df_loop <- cind_plot_df %>% 
  mutate(
    colour = case_when(
      dist == "exp" ~ "blue",
      dist == "gompertz" ~ "orange",
      dist == "llogis" ~ "deeppink",
      dist == "lnorm" ~ "purple",
      dist == "weibull" ~ "red"
    ),
    id = paste0(plot_id, " FOR ", cause)
  ) %>% 
  mutate(across(c(estimate, fit_est), ~ round((1 - .) * 100, 1)))

# Remove iterations per cause and trial causing plots to get compressed
df_loop_adj <- df_loop %>%
  left_join(thresholds, by = c("cause", "plot_id")) %>%
  group_by(plot_id, cause, iter) %>%
  filter(is.na(threshold) | all(fit_est <= threshold)) %>%
  ungroup()

# Check counts to make sure none were removed by accident
df_loop %>% 
  select(ctgov, cause, dist) %>%
  distinct() %>% 
  count(dist, cause) %>% 
  pivot_wider(names_from = "dist", values_from = "n")

df_loop_adj %>% 
  select(ctgov, cause, dist) %>%
  distinct() %>% 
  count(dist, cause) %>% 
  pivot_wider(names_from = "dist", values_from = "n")

# Remove unused objects to save memory
gc()
rm(imp_cind, cind_plot_df, df_loop)
gc()

# Get basic plots and blanks where a cause didn't occur
lst_plots <- get_plots()
gc()

# Figure 4 #####################################################################

# Prepare layout
layout_main <- df_loop_adj %>% 
  select(plot_id, cause) %>% 
  distinct() %>% 
  count(plot_id, cause) %>% 
  group_by(plot_id) %>% 
  reframe(n = length(unique(cause))) %>% 
  filter(n >= 5)

names_main <- paste0(unique(layout_main$plot_id), collapse = "|")

# Get plots
main_plots <- lst_plots[grepl(names_main, names(lst_plots))]

gt_fig4 <- dplyr::tibble(
  id = layout_main$plot_id,
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
    id = "Trial",
    ae_surv = "Adverse Event",
    loe_surv = "Lack of Efficacy",
    l2f_surv = "Lost to Follow-up",
    pi_surv = "PI/Sponsor Decision",
    pd_surv = "Protocol Violation",
    o_surv = "Other",
    vw_surv = "Voluntary Withdrawal"
  ) %>%
  tab_options(data_row.padding = px(1)) %>%
  tab_style(
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
      main_plots[grep("Adverse", names(main_plots))] %>%
        ggplot_image(height = px(30), aspect_ratio = 2)
    }
  ) %>%
  text_transform(  # Lack of Efficacy
    locations = cells_body(columns = loe_surv),
    fn = function(x) {
      main_plots[grep("Lack", names(main_plots))] %>%
        ggplot_image(height = px(30), aspect_ratio = 2)
    }
  ) %>%
  text_transform(  # Lost to Follow-up
    locations = cells_body(columns = l2f_surv),
    fn = function(x) {
      main_plots[grep("Lost", names(main_plots))] %>%
        ggplot_image(height = px(30), aspect_ratio = 2)
    }
  ) %>%
  text_transform(  # Other
    locations = cells_body(columns = o_surv),
    fn = function(x) {
      main_plots[grep("Other", names(main_plots))] %>%
        ggplot_image(height = px(30), aspect_ratio = 2)
    }
  ) %>%
  text_transform(  # PI/Sponsor Decision
    locations = cells_body(columns = pi_surv),
    fn = function(x) {
      main_plots[grep("PI", names(main_plots))] %>%
        ggplot_image(height = px(30), aspect_ratio = 2)
    }
  ) %>%
  text_transform(  # Protocol Deviation
    locations = cells_body(columns = pd_surv),
    fn = function(x) {
      main_plots[grep("Protocol", names(main_plots))] %>%
        ggplot_image(height = px(30), aspect_ratio = 2)
    }
  ) %>%
  text_transform(  # Voluntary Withdrawal
    locations = cells_body(columns = vw_surv),
    fn = function(x) {
      main_plots[grep("Voluntary", names(main_plots))] %>%
        ggplot_image(height = px(30), aspect_ratio = 2)
    }
  )

gtsave(gt_fig4, "categorical_attrition/plots/gt/main.html")

# Supp Figure 4a ###############################################################

layout_suppl <- df_loop_adj %>% 
  select(plot_id, cause) %>% 
  distinct() %>% 
  count(plot_id, cause) %>% 
  group_by(plot_id) %>% 
  reframe(n = length(unique(cause))) %>% 
  filter(n < 5)

layout_s4a <- layout_suppl[1:14,]
names_s4a <- paste0(unique(layout_s4a$plot_id), collapse = "|")
s4a_plots <- lst_plots[grepl(names_s4a, names(lst_plots))]

gt_s4a <- dplyr::tibble(
  id = layout_s4a$plot_id,
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
    id = "Trial",
    ae_surv = "Adverse Event",
    loe_surv = "Lack of Efficacy",
    l2f_surv = "Lost to Follow-up",
    pi_surv = "PI/Sponsor Decision",
    pd_surv = "Protocol Violation",
    o_surv = "Other",
    vw_surv = "Voluntary Withdrawal"
  ) %>%
  tab_options(data_row.padding = px(1)) %>%
  tab_style(
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
      s4a_plots[grep("Adverse", names(s4a_plots))] %>%
        ggplot_image(height = px(30), aspect_ratio = 2)
    }
  ) %>%
  text_transform(  # Lack of Efficacy
    locations = cells_body(columns = loe_surv),
    fn = function(x) {
      s4a_plots[grep("Lack", names(s4a_plots))] %>%
        ggplot_image(height = px(30), aspect_ratio = 2)
    }
  ) %>%
  text_transform(  # Lost to Follow-up
    locations = cells_body(columns = l2f_surv),
    fn = function(x) {
      s4a_plots[grep("Lost", names(s4a_plots))] %>%
        ggplot_image(height = px(30), aspect_ratio = 2)
    }
  ) %>%
  text_transform(  # Other
    locations = cells_body(columns = o_surv),
    fn = function(x) {
      s4a_plots[grep("Other", names(s4a_plots))] %>%
        ggplot_image(height = px(30), aspect_ratio = 2)
    }
  ) %>%
  text_transform(  # PI/Sponsor Decision
    locations = cells_body(columns = pi_surv),
    fn = function(x) {
      s4a_plots[grep("PI", names(s4a_plots))] %>%
        ggplot_image(height = px(30), aspect_ratio = 2)
    }
  ) %>%
  text_transform(  # Protocol Deviation
    locations = cells_body(columns = pd_surv),
    fn = function(x) {
      s4a_plots[grep("Protocol", names(s4a_plots))] %>%
        ggplot_image(height = px(30), aspect_ratio = 2)
    }
  ) %>%
  text_transform(  # Voluntary Withdrawal
    locations = cells_body(columns = vw_surv),
    fn = function(x) {
      s4a_plots[grep("Voluntary", names(s4a_plots))] %>%
        ggplot_image(height = px(30), aspect_ratio = 2)
    }
  )

gtsave(gt_s4a, "categorical_attrition/plots/gt/supp4a.html")

# Supp Figure 4b ###############################################################

layout_suppl <- df_loop_adj %>% 
  select(plot_id, cause) %>% 
  distinct() %>% 
  count(plot_id, cause) %>% 
  group_by(plot_id) %>% 
  reframe(n = length(unique(cause))) %>% 
  filter(n < 5)

layout_s4b <- layout_suppl[15:28,]
names_s4b <- paste0(unique(layout_s4b$plot_id), collapse = "|")
s4b_plots <- lst_plots[grepl(names_s4b, names(lst_plots))]

gt_s4b <- dplyr::tibble(
  id = layout_s4b$plot_id,
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
    id = "Trial",
    ae_surv = "Adverse Event",
    loe_surv = "Lack of Efficacy",
    l2f_surv = "Lost to Follow-up",
    pi_surv = "PI/Sponsor Decision",
    pd_surv = "Protocol Violation",
    o_surv = "Other",
    vw_surv = "Voluntary Withdrawal"
  ) %>%
  tab_options(data_row.padding = px(1)) %>%
  tab_style(
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
      s4b_plots[grep("Adverse", names(s4b_plots))] %>%
        ggplot_image(height = px(30), aspect_ratio = 2)
    }
  ) %>%
  text_transform(  # Lack of Efficacy
    locations = cells_body(columns = loe_surv),
    fn = function(x) {
      s4b_plots[grep("Lack", names(s4b_plots))] %>%
        ggplot_image(height = px(30), aspect_ratio = 2)
    }
  ) %>%
  text_transform(  # Lost to Follow-up
    locations = cells_body(columns = l2f_surv),
    fn = function(x) {
      s4b_plots[grep("Lost", names(s4b_plots))] %>%
        ggplot_image(height = px(30), aspect_ratio = 2)
    }
  ) %>%
  text_transform(  # Other
    locations = cells_body(columns = o_surv),
    fn = function(x) {
      s4b_plots[grep("Other", names(s4b_plots))] %>%
        ggplot_image(height = px(30), aspect_ratio = 2)
    }
  ) %>%
  text_transform(  # PI/Sponsor Decision
    locations = cells_body(columns = pi_surv),
    fn = function(x) {
      s4b_plots[grep("PI", names(s4b_plots))] %>%
        ggplot_image(height = px(30), aspect_ratio = 2)
    }
  ) %>%
  text_transform(  # Protocol Deviation
    locations = cells_body(columns = pd_surv),
    fn = function(x) {
      s4b_plots[grep("Protocol", names(s4b_plots))] %>%
        ggplot_image(height = px(30), aspect_ratio = 2)
    }
  ) %>%
  text_transform(  # Voluntary Withdrawal
    locations = cells_body(columns = vw_surv),
    fn = function(x) {
      s4b_plots[grep("Voluntary", names(s4b_plots))] %>%
        ggplot_image(height = px(30), aspect_ratio = 2)
    }
  )

gtsave(gt_s4b, "categorical_attrition/plots/gt/supp4b.html")

# Supp Figure 4c ###############################################################

layout_suppl <- df_loop_adj %>% 
  select(plot_id, cause) %>% 
  distinct() %>% 
  count(plot_id, cause) %>% 
  group_by(plot_id) %>% 
  reframe(n = length(unique(cause))) %>% 
  filter(n < 5)

layout_s4c <- layout_suppl[29:42,]
names_s4c <- paste0(unique(layout_s4c$plot_id), collapse = "|")
s4c_plots <- lst_plots[grepl(names_s4c, names(lst_plots))]

gt_s4c <- dplyr::tibble(
  id = layout_s4c$plot_id,
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
    id = "Trial",
    ae_surv = "Adverse Event",
    loe_surv = "Lack of Efficacy",
    l2f_surv = "Lost to Follow-up",
    pi_surv = "PI/Sponsor Decision",
    pd_surv = "Protocol Violation",
    o_surv = "Other",
    vw_surv = "Voluntary Withdrawal"
  ) %>%
  tab_options(data_row.padding = px(1)) %>%
  tab_style(
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
      s4c_plots[grep("Adverse", names(s4c_plots))] %>%
        ggplot_image(height = px(30), aspect_ratio = 2)
    }
  ) %>%
  text_transform(  # Lack of Efficacy
    locations = cells_body(columns = loe_surv),
    fn = function(x) {
      s4c_plots[grep("Lack", names(s4c_plots))] %>%
        ggplot_image(height = px(30), aspect_ratio = 2)
    }
  ) %>%
  text_transform(  # Lost to Follow-up
    locations = cells_body(columns = l2f_surv),
    fn = function(x) {
      s4c_plots[grep("Lost", names(s4c_plots))] %>%
        ggplot_image(height = px(30), aspect_ratio = 2)
    }
  ) %>%
  text_transform(  # Other
    locations = cells_body(columns = o_surv),
    fn = function(x) {
      s4c_plots[grep("Other", names(s4c_plots))] %>%
        ggplot_image(height = px(30), aspect_ratio = 2)
    }
  ) %>%
  text_transform(  # PI/Sponsor Decision
    locations = cells_body(columns = pi_surv),
    fn = function(x) {
      s4c_plots[grep("PI", names(s4c_plots))] %>%
        ggplot_image(height = px(30), aspect_ratio = 2)
    }
  ) %>%
  text_transform(  # Protocol Deviation
    locations = cells_body(columns = pd_surv),
    fn = function(x) {
      s4c_plots[grep("Protocol", names(s4c_plots))] %>%
        ggplot_image(height = px(30), aspect_ratio = 2)
    }
  ) %>%
  text_transform(  # Voluntary Withdrawal
    locations = cells_body(columns = vw_surv),
    fn = function(x) {
      s4c_plots[grep("Voluntary", names(s4c_plots))] %>%
        ggplot_image(height = px(30), aspect_ratio = 2)
    }
  )

gtsave(gt_s4c, "categorical_attrition/plots/gt/supp4c.html")

# Supp Figure 4d ###############################################################

layout_suppl <- df_loop_adj %>% 
  select(plot_id, cause) %>% 
  distinct() %>% 
  count(plot_id, cause) %>% 
  group_by(plot_id) %>% 
  reframe(n = length(unique(cause))) %>% 
  filter(n < 5)

layout_s4d <- layout_suppl[43:56,]
names_s4d <- paste0(unique(layout_s4d$plot_id), collapse = "|")
s4d_plots <- lst_plots[grepl(names_s4d, names(lst_plots))]

gt_s4d <- dplyr::tibble(
  id = layout_s4d$plot_id,
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
    id = "Trial",
    ae_surv = "Adverse Event",
    loe_surv = "Lack of Efficacy",
    l2f_surv = "Lost to Follow-up",
    pi_surv = "PI/Sponsor Decision",
    pd_surv = "Protocol Violation",
    o_surv = "Other",
    vw_surv = "Voluntary Withdrawal"
  ) %>%
  tab_options(data_row.padding = px(1)) %>%
  tab_style(
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
      s4d_plots[grep("Adverse", names(s4d_plots))] %>%
        ggplot_image(height = px(30), aspect_ratio = 2)
    }
  ) %>%
  text_transform(  # Lack of Efficacy
    locations = cells_body(columns = loe_surv),
    fn = function(x) {
      s4d_plots[grep("Lack", names(s4d_plots))] %>%
        ggplot_image(height = px(30), aspect_ratio = 2)
    }
  ) %>%
  text_transform(  # Lost to Follow-up
    locations = cells_body(columns = l2f_surv),
    fn = function(x) {
      s4d_plots[grep("Lost", names(s4d_plots))] %>%
        ggplot_image(height = px(30), aspect_ratio = 2)
    }
  ) %>%
  text_transform(  # Other
    locations = cells_body(columns = o_surv),
    fn = function(x) {
      s4d_plots[grep("Other", names(s4d_plots))] %>%
        ggplot_image(height = px(30), aspect_ratio = 2)
    }
  ) %>%
  text_transform(  # PI/Sponsor Decision
    locations = cells_body(columns = pi_surv),
    fn = function(x) {
      s4d_plots[grep("PI", names(s4d_plots))] %>%
        ggplot_image(height = px(30), aspect_ratio = 2)
    }
  ) %>%
  text_transform(  # Protocol Deviation
    locations = cells_body(columns = pd_surv),
    fn = function(x) {
      s4d_plots[grep("Protocol", names(s4d_plots))] %>%
        ggplot_image(height = px(30), aspect_ratio = 2)
    }
  ) %>%
  text_transform(  # Voluntary Withdrawal
    locations = cells_body(columns = vw_surv),
    fn = function(x) {
      s4d_plots[grep("Voluntary", names(s4d_plots))] %>%
        ggplot_image(height = px(30), aspect_ratio = 2)
    }
  )

gtsave(gt_s4d, "categorical_attrition/plots/gt/supp4d.html")

# Supp Figure 4e ###############################################################

layout_suppl <- df_loop_adj %>% 
  select(plot_id, cause) %>% 
  distinct() %>% 
  count(plot_id, cause) %>% 
  group_by(plot_id) %>% 
  reframe(n = length(unique(cause))) %>% 
  filter(n < 5)

layout_s4e <- layout_suppl[57:71,]
names_s4e <- paste0(unique(layout_s4e$plot_id), collapse = "|")
s4e_plots <- lst_plots[grepl(names_s4e, names(lst_plots))]

gt_s4e <- dplyr::tibble(
  id = layout_s4e$plot_id,
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
    id = "Trial",
    ae_surv = "Adverse Event",
    loe_surv = "Lack of Efficacy",
    l2f_surv = "Lost to Follow-up",
    pi_surv = "PI/Sponsor Decision",
    pd_surv = "Protocol Violation",
    o_surv = "Other",
    vw_surv = "Voluntary Withdrawal"
  ) %>%
  tab_options(data_row.padding = px(1)) %>%
  tab_style(
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
      s4e_plots[grep("Adverse", names(s4e_plots))] %>%
        ggplot_image(height = px(30), aspect_ratio = 2)
    }
  ) %>%
  text_transform(  # Lack of Efficacy
    locations = cells_body(columns = loe_surv),
    fn = function(x) {
      s4e_plots[grep("Lack", names(s4e_plots))] %>%
        ggplot_image(height = px(30), aspect_ratio = 2)
    }
  ) %>%
  text_transform(  # Lost to Follow-up
    locations = cells_body(columns = l2f_surv),
    fn = function(x) {
      s4e_plots[grep("Lost", names(s4e_plots))] %>%
        ggplot_image(height = px(30), aspect_ratio = 2)
    }
  ) %>%
  text_transform(  # Other
    locations = cells_body(columns = o_surv),
    fn = function(x) {
      s4e_plots[grep("Other", names(s4e_plots))] %>%
        ggplot_image(height = px(30), aspect_ratio = 2)
    }
  ) %>%
  text_transform(  # PI/Sponsor Decision
    locations = cells_body(columns = pi_surv),
    fn = function(x) {
      s4e_plots[grep("PI", names(s4e_plots))] %>%
        ggplot_image(height = px(30), aspect_ratio = 2)
    }
  ) %>%
  text_transform(  # Protocol Deviation
    locations = cells_body(columns = pd_surv),
    fn = function(x) {
      s4e_plots[grep("Protocol", names(s4e_plots))] %>%
        ggplot_image(height = px(30), aspect_ratio = 2)
    }
  ) %>%
  text_transform(  # Voluntary Withdrawal
    locations = cells_body(columns = vw_surv),
    fn = function(x) {
      s4e_plots[grep("Voluntary", names(s4e_plots))] %>%
        ggplot_image(height = px(30), aspect_ratio = 2)
    }
  )

gtsave(gt_s4e, "categorical_attrition/plots/gt/supp4e.html")

webshot(
  "categorical_attrition/plots/gt/supp4e.html", 
  file = "categorical_attrition/plots/gt/supp4e.png", 
  vwidth = 2244, 
  vheight = 1890, 
  zoom = 2
)













