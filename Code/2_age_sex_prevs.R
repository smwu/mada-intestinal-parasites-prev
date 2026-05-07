#===============================================================================
# Intestinal Parasites Age-Sex Prevalence Estimates
# Author: SM Wu
# Date created: 2025/01/28
# Date updated: 2025/11/24
# Purpose: Create exploratory figures and parasite prevalence estimates 
# STEPS: 
# (1) Read in data
# (2) Get age-sex prevalence estimate
# (3) Plot age-sex-specific prevalences 
# (4) Table of regression coefficients
# (5) Sanity check with manual calculations
# 
# Inputs:
#   Datasets:
#       "Model_Outputs/all_fec_bin_models_tp1.RData": Model outputs for all 
#          parasites using data from the first available timepoint
#  
# Outputs (number in parentheses indicates step in which it was generated):
#   Final outputs: 
#       (2) "Model_Outputs/prev_age_sex_all.csv": Dataframe of age-sex prevalences
#       (2) "Model_Outputs/prev_age_sex_table.csv": Table of age-sex prevalences
#       (3) "Model_Outputs/Figures/age_sex_heatmap_facetted.png": Plot of region prevalences for all parasites, sharing a color scale
#       (3) "Model_Outputs/Figures/age_sex_heatmap_sep.png": Plot of region prevalences for all parasites, with separate color scales per parasite


#============== (1) Read in data ===============================================
# Clear memory
rm(list = ls())

# Load libraries
library(tidyverse)        # for tidy data routines
library(magrittr)         # for efficient piping
library(ggplot2)          # for plotting
library(ggbeeswarm)       # for beeswarm plots
library(brms)             # for Bayesian hierarchical modeling
library(tidybayes)        # for Bayesian modeling helper functions
library(grid)             # for displaying tables
library(gridExtra)        # for displaying tables
library(writexl)          # for writing excel tables
library(broom)            # tidy data frames
library(gtsummary)        # create tables
library(flextable)        # tables in Word doc
library(officer)          # create Word docs
library(ggpubr)

### Specify directories (change to your local paths)
wd <- "~/Documents/GitHub/mada-intestinal-parasites-prev/"  # Working directory
setwd(wd)
data_dir <- "Cleaned_Data/"  # Directory with data 
code_dir <- "Code/"  # Directory with code
res_dir <- "Model_Outputs/"  # Directory to store results

# Read in prevalence models
fec_all <- load(paste0(wd, res_dir, "all_fec_bin_models_tp1.Rdata"))


#======== (2) Get age-sex prevalence estimates =================================


# Function to obtain age-sex specific prevalence estimates
get_age_sex_prev <- function(model, digits = 3) {
  # Define sex and age_cat levels
  new_data <- expand_grid(sex = c("M", "F"),
                          age_cat = c("[0,2)", "[2,5)", "[5,12)", "[12,20)", 
                                      "[20,50)", "[50,Inf)"))
  
  # Generate grand mean for age-sex groups (ignoring variability from random effects)
  grand_mean_samps <- model %>%
    epred_draws(newdata = new_data, allow_new_levels = TRUE, re_formula = NA) %>%
    rename(prevalence = .epred)
  
  # Create plot of prevalence distributions
  prevalence_plot <- grand_mean_samps %>% 
    ggplot(aes(x = prevalence, fill = sex)) +
    facet_grid(sex ~ age_cat) + 
    stat_halfeye() +
    labs(x = "Prevalence", y = NULL,
         fill = "Sex") +
    theme_bw() +
    theme(legend.position = "right")
  
  # Summarize results by age-sex category
  prevalence_summary <- grand_mean_samps %>% 
    select(sex, age_cat, prevalence) %>%
    group_by(sex, age_cat) %>%
    # pivot_longer(cols = everything(), names_to = "group", values_to = "prevalence") %>%
    # group_by(group) %>%
    summarise(
      mean = format(round(mean(prevalence), digits = digits), 
                    digits = digits, nsmall = digits),
      lower = format(round(quantile(prevalence, 0.025), digits = digits), 
                     digits = digits, nsmall = digits),
      upper = format(round(quantile(prevalence, 0.975), digits = digits), 
                     digits = digits, nsmall = digits)
    )
  
  # Reorder rows
  prevalence_summary <- 
    as.data.frame(
      prevalence_summary[c(7, 9, 11, 8, 10, 12, 1, 3, 5, 2, 4, 6), ])
  
  # Convert prevalences to numeric
  prevalence_summary <- prevalence_summary %>%
    mutate(mean = as.numeric(mean),
           lower = as.numeric(lower),
           upper = as.numeric(upper))
  
  # Print results
  return(list(prevalence_summary = prevalence_summary, 
              prevalence_plot = prevalence_plot))
}

# Get age-sex prevalence estimates for all parasites
# Use main models containing age_cat, sex fixed effects
### (1) Ascaris
prev_age_sex_ascaris <- get_age_sex_prev(model = ascaris_bin_main)
### (2) Hookworm
prev_age_sex_hook <- get_age_sex_prev(model = hook_bin_main)
### (3) Trichuris
prev_age_sex_trich <- get_age_sex_prev(model = trich_bin_main)
### (4) Strongyloides
prev_age_sex_strongyloides <- get_age_sex_prev(model = strongyloides_bin_main)
### (5) H.nana
prev_age_sex_h_nana <- get_age_sex_prev(model = h_nana_bin_main)
### (6) S.mansoni
prev_age_sex_s_mansoni <- get_age_sex_prev(model = s_mansoni_bin_main)
### (7) Helms
prev_age_sex_helms <- get_age_sex_prev(model = helms_bin_main)
### (8) E-coli
prev_age_sex_e_coli <- get_age_sex_prev(model = e_coli_bin_main)



# # Get age-sex prevalence estimates for all parasites
# # Use interaction models containing age_cat, sex, age_cat:sex fixed effects
# ### (1) Ascaris
# prev_age_sex_ascaris <- get_age_sex_prev(model = ascaris_bin_int)
# ### (2) Hookworm
# prev_age_sex_hook <- get_age_sex_prev(model = hook_bin_int)
# ### (3) Trichuris
# prev_age_sex_trich <- get_age_sex_prev(model = trich_bin_int)
# ### (4) Strongyloides
# prev_age_sex_strongyloides <- get_age_sex_prev(model = strongyloides_bin_int)
# ### (5) H.nana
# prev_age_sex_h_nana <- get_age_sex_prev(model = h_nana_bin_int)
# ### (6) S.mansoni
# prev_age_sex_s_mansoni <- get_age_sex_prev(model = s_mansoni_bin_int)
# ### (7) Helms
# prev_age_sex_helms <- get_age_sex_prev(model = helms_bin_int)
# ### (8) E-coli
# prev_age_sex_e_coli <- get_age_sex_prev(model = e_coli_bin_int)

# Append together all parasite prevalence estimates
prev_age_sex_all <- rbind(prev_age_sex_ascaris$prevalence_summary, 
                  prev_age_sex_hook$prevalence_summary, 
                  prev_age_sex_trich$prevalence_summary, 
                  prev_age_sex_strongyloides$prevalence_summary, 
                  prev_age_sex_h_nana$prevalence_summary, 
                  prev_age_sex_s_mansoni$prevalence_summary, 
                  prev_age_sex_helms$prevalence_summary, 
                  prev_age_sex_e_coli$prevalence_summary)
prev_age_sex_all$parasite <- rep(c("Ascaris", "Hookworm", "Trichuris", "Strongyloides", 
                                   "H. nana", "S. mansoni", "Helms", "E. coli"), 
                                 each = 12)

# # Display table
# grid::grid.newpage()
# gridExtra::grid.table(prev_age_sex_all)
# kable(prev_age_sex_all, format = "latex", booktabs = TRUE)

# # Save dataframe of age-sex prevalences
# write.csv(prev_age_sex_all, file = paste0(wd, res_dir, "Tables/prev_age_sex_all.csv"),
#           row.names = FALSE)

# Function to tidy up the prevalences to display in a table
tidy_prev <- function(prev_summ) {
  tidy_summ <- paste0(format(prev_summ$mean * 100, 1), " (", 
                      format(prev_summ$lower * 100, 1), ", ",
                      format(prev_summ$upper * 100, 1), ")")
  return(tidy_summ)
}

# Create nicer table
prev_age_sex_table <- data.frame(
  Sex = c(rep("Male", 6), rep("Female", 6)),
  Age = rep(c("<2", "2-4", "5-11", "12-19", "20-49", ">=50"), 2))
prev_age_sex_table$`A. lumbricoides` <- 
  tidy_prev(prev_age_sex_ascaris$prevalence_summary)
prev_age_sex_table$`T. trichiura` <- 
  tidy_prev(prev_age_sex_trich$prevalence_summary)
prev_age_sex_table$`Hookworm` <- 
  tidy_prev(prev_age_sex_hook$prevalence_summary)
prev_age_sex_table$`Strongyloides` <- 
  tidy_prev(prev_age_sex_strongyloides$prevalence_summary)
prev_age_sex_table$`H. nana` <- 
  tidy_prev(prev_age_sex_h_nana$prevalence_summary)
prev_age_sex_table$`S. mansoni` <- 
  tidy_prev(prev_age_sex_s_mansoni$prevalence_summary)
prev_age_sex_table$`Other Helminths` <- 
  tidy_prev(prev_age_sex_helms$prevalence_summary)
prev_age_sex_table$`E. coli` <- 
  tidy_prev(prev_age_sex_e_coli$prevalence_summary)

# # Save nice table of age and sex prevalences
# write_xlsx(prev_age_sex_table, 
#           path = paste0(wd, res_dir, "Tables/prev_age_sex_table.xlsx"))

#============== (3) Plot age-sex-specific prevalences ==========================

# Prepare plotting data
prev_age_sex_all_map <- prev_age_sex_all %>%
  rename(Prevalence = mean,
         Age = age_cat,
         Sex = sex,
         Parasite = parasite) %>%
  mutate(Sex = factor(Sex, levels = c("M", "F")),
         Age = factor(Age, levels = c("[0,2)", "[2,5)", "[5,12)", "[12,20)", 
                                      "[20,50)", "[50,Inf)"),
                      labels = c("<2", "2-4", "5-11", "12-19", "20-49", ">=50")),
         Parasite = factor(Parasite, 
                           levels = c("Ascaris", "Trichuris", "Hookworm", 
                                                "Strongyloides", "H. nana", 
                                                "S. mansoni", "Helms", "E. coli"),
                           labels = c("A.lumbricoides", "T.trichiura", "Hookworm", 
                                      "Strongyloides", "H.nana", "S.mansoni", 
                                      "Helminths", "E.coli")))
parasites <- c("A.lumbricoides", "T.trichiura", "Hookworm", 
               "Strongyloides", "H.nana", "S.mansoni", 
               "Helminths", "E.coli")

# Plot of age-sex prevalences
### Creating heatmap of region-specific prevalences for all parasites within a 
### single facetted plot
prev_age_sex_all_map %>% 
  mutate(Sex = factor(Sex, levels = c("M", "F"), labels = c("Males", "Females"))) %>%
  ggplot(aes(x = Parasite, y = fct_rev(Age), fill = (Prevalence * 100))) + 
  facet_grid(Sex ~ Parasite, scales = "free") + 
  geom_tile(color = "white",
            lwd = 1.5,
            linetype = 1) +
  scale_fill_gradientn(
    colours = c("#FFFEF5", "#FFF6D5", "#FEECC2", "#FED7A2", "#FDBB84", "#FC8D59", "#D7301F", "#8F477D", "#612DA4"),
    values  = scales::rescale(c(0, 0.2, 0.5, 1, 5, 20, 35, 55, 65)),
    limits  = c(0, 65),
    oob     = scales::squish,
    na.value = "white"
  ) + 
  geom_text(aes(label = round(Prevalence * 100, 1)), 
            color = "black", size = 4) +
  theme_bw() +
  theme(panel.spacing = unit(0.4, "cm"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 11),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.key.width = unit(2, "cm"),
        strip.background = element_rect(fill = "snow2"),
        strip.text = element_text(face = "bold", size = 11)) + 
  guides(fill = guide_colourbar(title = "Prevalence (%)")) +
  labs(y = "Age")
# # Save plot
# ggsave(filename = paste0(wd, res_dir, "Figures/age_sex_heatmap_facetted.png"),
#        width = 11.5, height = 7.5, units = "in")

# # Export as EPS to convert to TIFF with NAAS tool
# ggsave(filename = paste0(wd, res_dir, "Figures/Fig3.eps"),
#        width = 11.5, height = 7.5, units = "in")



### Creating heatmap of age-sex-specific prevalences for all parasites, each with 
### its own legend and axes
get_age_sex_heatmap <- function(prev_data, parasite_index) {
  heatmap <- prev_age_sex_all_map %>% 
    filter(Parasite == parasites[parasite_index]) %>% 
    mutate(Sex = factor(Sex, levels = c("M", "F"), labels = c("Males", "Females"))) %>%
    ggplot(aes(x = Parasite, y = fct_rev(Age), fill = (Prevalence * 100))) + 
    coord_fixed() + 
    facet_grid(Sex ~ .) + 
    geom_tile(color = "white",
              lwd = 1.5,
              linetype = 1) +
    scale_fill_gradient(low = "#fef0d9", high = "#d7301f", na.value = "white",
                        labels = function(x) format(x, digits = 1)) +
    geom_text(aes(label = format(round(Prevalence * 100, 1), 
                                 digits = 1, nsmall = 1)), 
              color = "black", size = 4) +
    theme_bw() +
    theme(panel.spacing = unit(0.1, "cm")) + 
    theme(plot.margin = unit(c(0, 0.5, 0, 0.2), 'lines')) + 
    theme(strip.background = element_rect(fill = "snow2")) +
    theme(axis.title = element_blank(),
          legend.position = "bottom",
          legend.direction = "horizontal") + 
    guides(fill = guide_colourbar(title = "Prev\n (%)")) +
    # guides(fill = guide_colourbar(title = "")) +
    theme(
      # strip.text = element_text(face = "bold"),
      axis.text.x = element_text(size = 10, face = "bold"),  # Bold parasite names
      axis.text.y = element_text(size = 9),
      strip.text = element_text(size = 10),
      legend.title = element_text(size = 8),        # Smaller legend title
      legend.text = element_text(size = 7),         # Smaller legend text
      legend.key.size = unit(0.4, "cm"),            # Smaller legend keys
      legend.key.width = unit(0.4, "cm")           # Smaller key width (for better spacing)
    ) + 
    labs(y = "Age")
  
  return(heatmap)
}

# Get plots for all parasites and append together
p1 <- get_age_sex_heatmap(prev_data = prev_reg_all_map, parasite_index = 1)
p2 <- get_age_sex_heatmap(prev_data = prev_reg_all_map, parasite_index = 2)
p3 <- get_age_sex_heatmap(prev_data = prev_reg_all_map, parasite_index = 3)
p4 <- get_age_sex_heatmap(prev_data = prev_reg_all_map, parasite_index = 4)
p5 <- get_age_sex_heatmap(prev_data = prev_reg_all_map, parasite_index = 5)
p6 <- get_age_sex_heatmap(prev_data = prev_reg_all_map, parasite_index = 6)
p7 <- get_age_sex_heatmap(prev_data = prev_reg_all_map, parasite_index = 7)
p8 <- get_age_sex_heatmap(prev_data = prev_reg_all_map, parasite_index = 8)

p_all <- ggarrange(plotlist = list(p1, p2, p3, p4, p5, p6, p7, p8), 
                   ncol = 8, nrow = 1)
p_all <- annotate_figure(p_all, left = textGrob("Age", rot = 90, vjust = 1, 
                                       gp = gpar(cex = 1)))
# # Save plot of age-sex prevalences
# ggsave(filename = paste0(wd, res_dir, "Figures/age_sex_heatmap_sep.png"),
#        plot = p_all, width = 11, height = 7, units = "in")



### Create plot of age-sex prevalences and their CIs for all parasites
# Males
p1 <- prev_age_sex_all_map %>%
  filter(Sex == "M") %>%
  ggplot(aes(x = Prevalence, y = fct_rev(Age), xmin = lower, xmax = upper)) + 
  geom_errorbar(col = "darkgrey") + 
  geom_point(fill = "lightblue2", pch = 21, size = 3) + 
  theme_bw() +
  facet_wrap(~ Parasite, nrow = 1, ncol = 8, scales = "fixed") + 
  theme(panel.spacing = unit(0.3, "cm"),
        strip.text = element_text(size = 11)) + 
  theme(strip.background = element_rect(fill = "lightblue2")) +
  ylab("Age") + xlab("Prevalence for Males") + xlim(c(0, 1))
# Females
p2 <- prev_age_sex_all_map %>%
  filter(Sex == "F") %>%
  ggplot(aes(x = Prevalence, y = fct_rev(Age), xmin = lower, xmax = upper)) + 
  geom_errorbar(col = "darkgrey") + 
  geom_point(fill = "lightyellow", pch = 21, size = 3) + 
  theme_bw() +
  facet_wrap(~ Parasite, nrow = 1, ncol = 8, scales = "fixed") + 
  theme(panel.spacing = unit(0.3, "cm"),
        strip.text = element_text(size = 11)) + 
  theme(strip.background = element_rect(fill = "lightyellow")) +
  ylab("Age") + xlab("Prevalence for Females") + xlim(c(0, 1))
p_age_sex <- gridExtra::grid.arrange(arrangeGrob(p1, p2, nrow = 2))
# # Save plot
# ggsave(filename = paste0(wd, res_dir, "Figures/prev_age_sex_ci.png"),
#        plot = p_age_sex, width = 14, height = 7, units = "in")


#=========================
### Repeat for interaction models containing age_cat, sex, age_cat:sex fixed effects
### (1) Ascaris
prev_age_sex_ascaris <- get_age_sex_prev(model = ascaris_bin_int)
### (2) Hookworm
prev_age_sex_hook <- get_age_sex_prev(model = hook_bin_int)
### (3) Trichuris
prev_age_sex_trich <- get_age_sex_prev(model = trich_bin_int)
### (4) Strongyloides
prev_age_sex_strongyloides <- get_age_sex_prev(model = strongyloides_bin_int)
### (5) H.nana
prev_age_sex_h_nana <- get_age_sex_prev(model = h_nana_bin_int)
### (6) S.mansoni
prev_age_sex_s_mansoni <- get_age_sex_prev(model = s_mansoni_bin_int)
### (7) Helms
prev_age_sex_helms <- get_age_sex_prev(model = helms_bin_int)
### (8) E-coli
prev_age_sex_e_coli <- get_age_sex_prev(model = e_coli_bin_int)

# Append together all parasite prevalence estimates
prev_age_sex_all <- rbind(prev_age_sex_ascaris$prevalence_summary, 
                          prev_age_sex_hook$prevalence_summary, 
                          prev_age_sex_trich$prevalence_summary, 
                          prev_age_sex_strongyloides$prevalence_summary, 
                          prev_age_sex_h_nana$prevalence_summary, 
                          prev_age_sex_s_mansoni$prevalence_summary, 
                          prev_age_sex_helms$prevalence_summary, 
                          prev_age_sex_e_coli$prevalence_summary)
prev_age_sex_all$parasite <- rep(c("Ascaris", "Hookworm", "Trichuris", "Strongyloides", 
                                   "H. nana", "S. mansoni", "Helms", "E. coli"), 
                                 each = 12)

# # Display table
# grid::grid.newpage()
# gridExtra::grid.table(prev_age_sex_all)
# kable(prev_age_sex_all, format = "latex", booktabs = TRUE)

# # Save dataframe of age-sex prevalences: interactions
# write.csv(prev_age_sex_all, file = paste0(wd, res_dir, "Tables/prev_age_sex_all_int.csv"),
#           row.names = FALSE)

# Create nicer table: interactions
prev_age_sex_table <- data.frame(
  Sex = c(rep("Male", 6), rep("Female", 6)),
  Age = rep(c("<2", "2-4", "5-11", "12-19", "20-49", ">=50"), 2))
prev_age_sex_table$`A. lumbricoides` <- 
  tidy_prev(prev_age_sex_ascaris$prevalence_summary)
prev_age_sex_table$`T. trichiura` <- 
  tidy_prev(prev_age_sex_trich$prevalence_summary)
prev_age_sex_table$`Hookworm` <- 
  tidy_prev(prev_age_sex_hook$prevalence_summary)
prev_age_sex_table$`Strongyloides` <- 
  tidy_prev(prev_age_sex_strongyloides$prevalence_summary)
prev_age_sex_table$`H. nana` <- 
  tidy_prev(prev_age_sex_h_nana$prevalence_summary)
prev_age_sex_table$`S. mansoni` <- 
  tidy_prev(prev_age_sex_s_mansoni$prevalence_summary)
prev_age_sex_table$`Other Helminths` <- 
  tidy_prev(prev_age_sex_helms$prevalence_summary)
prev_age_sex_table$`E. coli` <- 
  tidy_prev(prev_age_sex_e_coli$prevalence_summary)

# # Save nice table of age and sex prevalences: interactions
# write_xlsx(prev_age_sex_table, 
#           path = paste0(wd, res_dir, "Tables/prev_age_sex_table_int.xlsx"))


## Plot age-sex-specific prevalences 
# Prepare plotting data
prev_age_sex_all_map <- prev_age_sex_all %>%
  rename(Prevalence = mean,
         Age = age_cat,
         Sex = sex,
         Parasite = parasite) %>%
  mutate(Sex = factor(Sex, levels = c("M", "F")),
         Age = factor(Age, levels = c("[0,2)", "[2,5)", "[5,12)", "[12,20)", 
                                      "[20,50)", "[50,Inf)"),
                      labels = c("<2", "2-4", "5-11", "12-19", "20-49", ">=50")),
         Parasite = factor(Parasite, 
                           levels = c("Ascaris", "Trichuris", "Hookworm", 
                                      "Strongyloides", "H. nana", 
                                      "S. mansoni", "Helms", "E. coli"),
                           labels = c("A.lumbricoides", "T.trichiura", "Hookworm", 
                                      "Strongyloides", "H.nana", "S.mansoni", 
                                      "Helminths", "E.coli")))
parasites <- c("A.lumbricoides", "T.trichiura", "Hookworm", 
               "Strongyloides", "H.nana", "S.mansoni", 
               "Helminths", "E.coli")

# Plot of age-sex prevalences
### Creating heatmap of region-specific prevalences for all parasites within a 
### single facetted plot: interactions
prev_age_sex_all_map %>% 
  ggplot(aes(x = Parasite, y = fct_rev(Age), fill = (Prevalence * 100))) + 
  facet_grid(Sex ~ Parasite, scales = "free") + 
  geom_tile(color = "white",
            lwd = 1.5,
            linetype = 1) +
  scale_fill_gradient(low = "#fef0d9", high = "#d7301f", na.value = "white") +
  geom_text(aes(label = round(Prevalence * 100, 1)), 
            color = "black", size = 4) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "snow2")) +
  theme(panel.spacing = unit(0.3, "cm")) + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "bottom",
        legend.direction = "horizontal") + 
  guides(fill = guide_colourbar(title = "Prevalence (%)")) +
  labs(y = "Age")
# # Save plot
# ggsave(filename = paste0(wd, res_dir, "Figures/age_sex_heatmap_facetted_int.png"),
#        width = 8.5, height = 7.5, units = "in")


### Creating heatmap of age-sex-specific prevalences for all parasites, each with 
### its own legend and axes: interactions
get_age_sex_heatmap <- function(prev_data, parasite_index) {
  heatmap <- prev_age_sex_all_map %>% 
    filter(Parasite == parasites[parasite_index]) %>% 
    mutate(Sex = factor(Sex, levels = c("M", "F"), labels = c("Males", "Females"))) %>%
    ggplot(aes(x = Parasite, y = fct_rev(Age), fill = (Prevalence * 100))) + 
    coord_fixed() + 
    facet_grid(Sex ~ .) + 
    geom_tile(color = "white",
              lwd = 1.5,
              linetype = 1) +
    scale_fill_gradient(low = "#fef0d9", high = "#d7301f", na.value = "white",
                        labels = function(x) format(x, digits = 1)) +
    geom_text(aes(label = format(round(Prevalence * 100, 1), 
                                 digits = 1, nsmall = 1)), 
              color = "black", size = 4) +
    theme_bw() +
    theme(panel.spacing = unit(0.1, "cm")) + 
    theme(plot.margin = unit(c(0, 0.5, 0, 0.2), 'lines')) + 
    theme(strip.background = element_rect(fill = "snow2")) +
    theme(axis.title = element_blank(),
          legend.position = "bottom",
          legend.direction = "horizontal") + 
    guides(fill = guide_colourbar(title = "Prev\n (%)")) +
    # guides(fill = guide_colourbar(title = "")) +
    theme(
      # strip.text = element_text(face = "bold"),
      axis.text.x = element_text(size = 10, face = "bold"),  # Bold parasite names
      axis.text.y = element_text(size = 9),
      strip.text = element_text(size = 10),
      legend.title = element_text(size = 8),        # Smaller legend title
      legend.text = element_text(size = 7),         # Smaller legend text
      legend.key.size = unit(0.4, "cm"),            # Smaller legend keys
      legend.key.width = unit(0.4, "cm")           # Smaller key width (for better spacing)
    ) + 
    labs(y = "Age")
  
  return(heatmap)
}

# Get plots for all parasites and append together
p1 <- get_age_sex_heatmap(prev_data = prev_reg_all_map, parasite_index = 1)
p2 <- get_age_sex_heatmap(prev_data = prev_reg_all_map, parasite_index = 2)
p3 <- get_age_sex_heatmap(prev_data = prev_reg_all_map, parasite_index = 3)
p4 <- get_age_sex_heatmap(prev_data = prev_reg_all_map, parasite_index = 4)
p5 <- get_age_sex_heatmap(prev_data = prev_reg_all_map, parasite_index = 5)
p6 <- get_age_sex_heatmap(prev_data = prev_reg_all_map, parasite_index = 6)
p7 <- get_age_sex_heatmap(prev_data = prev_reg_all_map, parasite_index = 7)
p8 <- get_age_sex_heatmap(prev_data = prev_reg_all_map, parasite_index = 8)

p_all <- ggarrange(plotlist = list(p1, p2, p3, p4, p5, p6, p7, p8), 
                   ncol = 8, nrow = 1)
p_all <- annotate_figure(p_all, left = textGrob("Age", rot = 90, vjust = 1, 
                                                gp = gpar(cex = 1)))
# # Save plot of age-sex prevalences: interactions
# ggsave(filename = paste0(wd, res_dir, "Figures/age_sex_heatmap_sep_int.png"),
#        plot = p_all, width = 11, height = 7, units = "in")



### Create plot of age-sex prevalences and their CIs for all parasites: interactions
# Males
p1 <- prev_age_sex_all_map %>%
  filter(Sex == "M") %>%
  ggplot(aes(x = Prevalence, y = fct_rev(Age), xmin = lower, xmax = upper)) + 
  geom_errorbar(col = "darkgrey") + 
  geom_point(fill = "lightblue2", pch = 21, size = 3) + 
  theme_bw() +
  facet_wrap(~ Parasite, nrow = 1, ncol = 8, scales = "fixed") + 
  theme(panel.spacing = unit(0.3, "cm"),
        strip.text = element_text(size = 11)) + 
  theme(strip.background = element_rect(fill = "lightblue2")) +
  ylab("Age") + xlab("Prevalence for Males") + xlim(c(0, 1))
# Females
p2 <- prev_age_sex_all_map %>%
  filter(Sex == "F") %>%
  ggplot(aes(x = Prevalence, y = fct_rev(Age), xmin = lower, xmax = upper)) + 
  geom_errorbar(col = "darkgrey") + 
  geom_point(fill = "lightyellow", pch = 21, size = 3) + 
  theme_bw() +
  facet_wrap(~ Parasite, nrow = 1, ncol = 8, scales = "fixed") + 
  theme(panel.spacing = unit(0.3, "cm"),
        strip.text = element_text(size = 11)) + 
  theme(strip.background = element_rect(fill = "lightyellow")) +
  ylab("Age") + xlab("Prevalence for Females") + xlim(c(0, 1))
p_age_sex <- gridExtra::grid.arrange(arrangeGrob(p1, p2, nrow = 2))
# # Save plot
# ggsave(filename = paste0(wd, res_dir, "Figures/prev_age_sex_ci_int.png"),
#        plot = p_age_sex, width = 14, height = 7, units = "in")


#============== (4) Table of regression coefficients ===========================
models <- list(`Ascaris lumbricoides` = ascaris_bin_int, 
               `Trichurs trichiura` = trich_bin_int, 
               `Hookworm` = hook_bin_int, 
               `Strongyloides` = strongyloides_bin_int, 
               `Hymenolepis nana` = h_nana_bin_int, 
               `Schistosoma mansoni` = s_mansoni_bin_int, 
               `Other helminths` = helms_bin_int, 
               `Escherichia coli` = e_coli_bin_int)  

# Create a Word document
doc <- read_docx()

# Loop through models and add their summaries to the document
for (model_name in names(models)) {
  model <- models[[model_name]]
  
  tbl <- summary(model)$fixed
  tbl_nice <- tbl[, c(1,3,4)]
  tbl_nice <- exp(tbl_nice) # Exponentiate to OR scale
  rownames(tbl_nice) <- c("Intercept", "SexMale", "Age<2", "Age2-4", "Age5-11",
                          "Age12-19", "Age>=50", 
                          "SexMale:Age<2", "SexMale:Age2-4", 
                          "SexMale:Age5-11", "SexMale:Age12-19", 
                          "SexMale:Age>=50")
  tbl_nice <- format(round(tbl_nice, 2), 2)
  ft <- flextable(tbl_nice) %>%
    autofit() %>%
    theme_booktabs() %>%
    set_table_properties(width = 1, layout = "autofit")
  
  # Add model name as a heading
  doc <- doc %>%
    body_add_par(paste(model_name), style = "heading 1") %>%
    body_add_flextable(ft) %>%
    body_add_par("")  # Add spacing
}

# # Save the Word file
# print(doc, target = paste0(wd, res_dir, "Tables/model_summaries.docx"))


#============== (5) Sanity check with manual calculations ======================

# # Function to obtain age-sex specific prevalence estimates manually
# get_age_sex_prev_manual <- function(model, digits = 3) {
#   model <- ascaris_bin_main  
#   posterior_samples <- as_draws_df(model)
#   
#   # Initialize data frame of results 
#   prev_df <- as.data.frame(matrix(NA, nrow = (n_age*n_sex), ncol = 7))
#   colnames(prev_df) <- c("Sex", "Age", "Mean", "SD", "2.5%", "50%", "97.5%")
#   prev_df$Sex <- c(rep("F", n_age), rep("M", n_age))
#   prev_df$Age <- rep(age_groups, times = n_sex)
#   
#   ## Females
#   # Ages [0,2)
#   prev <- plogis(rowSums(posterior_samples[, c("b_Intercept", "b_age_cat02")]))
#   prev_df[1, -c(1:2)] <- c(mean = mean(prev), sd = sd(prev), 
#                            quantile(prev, probs = c(0.025, 0.5, 0.975)))
#   # Ages [2,5)
#   prev <- plogis(rowSums(posterior_samples[, c("b_Intercept", "b_age_cat25")]))
#   prev_df[2, -c(1:2)] <- c(mean = mean(prev), sd = sd(prev), 
#                            quantile(prev, probs = c(0.025, 0.5, 0.975)))
#   # Ages [5,12)
#   prev <- plogis(rowSums(posterior_samples[, c("b_Intercept", "b_age_cat512")]))
#   prev_df[3, -c(1:2)] <- c(mean = mean(prev), sd = sd(prev), 
#                            quantile(prev, probs = c(0.025, 0.5, 0.975)))
#   # Ages [12,20)
#   prev <- plogis(rowSums(posterior_samples[, c("b_Intercept", "b_age_cat1220")]))
#   prev_df[4, -c(1:2)] <- c(mean = mean(prev), sd = sd(prev), 
#                            quantile(prev, probs = c(0.025, 0.5, 0.975)))
#   # Ages [20,50)
#   prev <- plogis(rowSums(posterior_samples[, c("b_Intercept")]))
#   prev_df[5, -c(1:2)] <- c(mean = mean(prev), sd = sd(prev), 
#                            quantile(prev, probs = c(0.025, 0.5, 0.975)))
#   # Ages [50,Inf)
#   prev <- plogis(rowSums(posterior_samples[, c("b_Intercept", "b_age_cat50Inf")]))
#   prev_df[6, -c(1:2)] <- c(mean = mean(prev), sd = sd(prev), 
#                            quantile(prev, probs = c(0.025, 0.5, 0.975)))
#   ## Males
#   # Ages [0,2)
#   prev <- plogis(rowSums(posterior_samples[, c("b_Intercept", "b_sexM", "b_age_cat02")]))
#   prev_df[7, -c(1:2)] <- c(mean = mean(prev), sd = sd(prev), 
#                            quantile(prev, probs = c(0.025, 0.5, 0.975)))
#   # Ages [2,5)
#   prev <- plogis(rowSums(posterior_samples[, c("b_Intercept", "b_sexM", "b_age_cat25")]))
#   prev_df[8, -c(1:2)] <- c(mean = mean(prev), sd = sd(prev), 
#                            quantile(prev, probs = c(0.025, 0.5, 0.975)))
#   # Ages [5,12)
#   prev <- plogis(rowSums(posterior_samples[, c("b_Intercept", "b_sexM", "b_age_cat512")]))
#   prev_df[9, -c(1:2)] <- c(mean = mean(prev), sd = sd(prev), 
#                            quantile(prev, probs = c(0.025, 0.5, 0.975)))
#   # Ages [12,20)
#   prev <- plogis(rowSums(posterior_samples[, c("b_Intercept", "b_sexM", "b_age_cat1220")]))
#   prev_df[10, -c(1:2)] <- c(mean = mean(prev), sd = sd(prev), 
#                             quantile(prev, probs = c(0.025, 0.5, 0.975)))
#   # Ages [20,50)
#   prev <- plogis(rowSums(posterior_samples[, c("b_Intercept", "b_sexM")]))
#   prev_df[11, -c(1:2)] <- c(mean = mean(prev), sd = sd(prev), 
#                             quantile(prev, probs = c(0.025, 0.5, 0.975)))
#   # Ages [50,Inf)
#   prev <- plogis(rowSums(posterior_samples[, c("b_Intercept", "b_sexM", "b_age_cat50Inf")]))
#   prev_df[12, -c(1:2)] <- c(mean = mean(prev), sd = sd(prev), 
#                             quantile(prev, probs = c(0.025, 0.5, 0.975)))
#   
#   prev_df[, -c(1:2)] <- format(round(prev_df[, -c(1:2)], digits), 
#                                digits = digits, nsmall = digits)
#   
#   return(prev_df)
# }
# 
# 
# # Compare prevalences with manual calculations
# # Careful: ordering of rows differs
# prev_age_sex_ascaris <- get_age_sex_prev(model = ascaris_bin_main)
# prev_age_sex_ascaris_manual <- get_age_sex_prev_manual(model = ascaris_bin_main)
# 
# prev_age_sex_ascaris
# prev_age_sex_ascaris_manual
# 
# 
# # Compare with crude estimates
# sapply(c("[0,2)", "[2,5)", "[5,12)", "[12,20)", "[20,50)", "[50,Inf)"), 
#        function(x) round(mean(fec_key_bin$ascaris_bin[fec_key_bin$sex == "F" & 
#                                                         fec_key_bin$age_cat == x], 
#                               na.rm = TRUE), 3))
# sapply(c("[0,2)", "[2,5)", "[5,12)", "[12,20)", "[20,50)", "[50,Inf)"), 
#        function(x) round(mean(fec_key_bin$ascaris_bin[fec_key_bin$sex == "M" & 
#                                                         fec_key_bin$age_cat == x], 
#                               na.rm = TRUE), 3))
# 
