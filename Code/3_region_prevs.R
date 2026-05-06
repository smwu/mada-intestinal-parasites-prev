#===============================================================================
# Intestinal Parasites Region-Specific Prevalence Estimates
# Author: SM Wu
# Date created: 2025/02/05
# Date updated: 2026/12/16
# Purpose: Create region-specific parasite prevalence estimates and plots
# STEPS: 
# (1) Read in data
# (2) Calculate region-specific prevalences
# (3) Plot region-specific prevalences 
# (4) Sanity check with manual calculations
# 
# Inputs:
#   Datasets:
#       "Cleaned_Data/fec_key_bin_tp1.csv": Data with binary parasite outcome used for models, subsetted to tp1
#       "Model_Outputs/all_fec_bin_models_tp1.RData": Model outputs for all parasites using tp1 data
#  
# Outputs (number in parentheses indicates step in which it was generated):
#   Final outputs: 
#       (2) "Model_Outputs/prev_reg_all.csv": Dataframe of region prevalences
#       (3) "Model_Outputs/Figures/region_heatmap_facetted.png": Plot of region prevalences for all parasites, sharing a color scale
#       (3) "Model_Outputs/Figures/region_heatmap_sep.png": Plot of region prevalences for all parasites, with separate color scales per parasite
#       (3) "Model_Outputs/Figures/madagascar_map.png": Plot of Madagascar map with site locations
#       (3) "Model_Outputs/Figures/prev_region_ci.png": Plot or region prevalences and 95% credible intervals for all parasites


#============== (1) Read in data ===============================================
# Clear memory
rm(list = ls())

# Load libraries
library(tidyverse)        # for tidy data routines
library(ggplot2)          # for plotting
library(grid)             # for displaying tables
library(gridExtra)        # for displaying tables
library(ggpubr)           # for plotting
library(ggthemes)         # for plotting
library(reshape2)         # for plotting
library(scales)           # for plotting
library(patchwork)        # for multiple plots
library(cowplot)          # for shared legend
library(tidybayes)        # for Bayesian modeling helper functions
library(brms)
library(writexl)
library(readxl)
library(maps)             # for Madagascar map
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

### Specify directories (change to your local paths)
wd <- "~/Documents/GitHub/mada-intestinal-parasites-prev/"  # Working directory
setwd(wd)
data_dir <- "Cleaned_Data/"  # Directory with data 
code_dir <- "Code/"  # Directory with code
res_dir <- "Model_Outputs/"  # Directory to store results

# Read in dataset used
fec_key_bin <- read.csv(paste0(wd, data_dir, "fec_key_bin_tp1.csv"))

# Read in prevalence models
fec_all <- load(paste0(wd, res_dir, "all_fec_bin_models_tp1.Rdata"))


#=============== (2) Calculate region-specific prevalences =====================

# Function to obtain region specific prevalence estimates
get_reg_prev <- function(model, digits = 3) {
  # Define region levels
  new_data <- data.frame(region = factor(c(1, 2, 3, 4, 5)))
  
  # Generate regional samples
  region_samps <- model %>%
    epred_draws(newdata = new_data, allow_new_levels = TRUE, 
                re_formula = ~ (1 | region)) %>%
    rename(prevalence = .epred)
  
  # Create plot of prevalence distributions
  prevalence_plot <- region_samps %>% 
    ggplot(aes(x = prevalence, fill = region)) +
    facet_grid(. ~ region) + 
    stat_halfeye() +
    labs(x = "Prevalence", y = NULL,
         fill = "Region") +
    theme_bw() +
    theme(legend.position = "right")
  
  # Summarize results by age-sex category
  prevalence_summary <- region_samps %>% 
    select(region, prevalence) %>%
    group_by(region) %>%
    summarise(
      mean = format(round(mean(prevalence), digits = digits), 
                    digits = digits, nsmall = digits),
      lower = format(round(quantile(prevalence, 0.025), digits = digits), 
                     digits = digits, nsmall = digits),
      upper = format(round(quantile(prevalence, 0.975), digits = digits), 
                     digits = digits, nsmall = digits)
    )
  # Relabel regions
  prevalence_summary$region <- case_match(prevalence_summary$region, 
                                          "1" ~ "NE", 
                                          "2" ~ "SE", 
                                          "3" ~ "SW", 
                                          "4" ~ "WC", 
                                          "5" ~ "CP",
                                          .default = NA)
  # Convert prevalences to numeric
  prevalence_summary <- prevalence_summary %>%
    mutate(mean = as.numeric(mean),
           lower = as.numeric(lower),
           upper = as.numeric(upper))
  
  # Print results
  return(list(prevalence_summary = prevalence_summary, 
              prevalence_plot = prevalence_plot))
}


# Function to obtain prevalence estimates from brms model when only 1 region has
# data available. We don't want to predict for other regions because we don't 
# have enough data to estimate a full distribution of region-level effects
get_prev_1_reg <- function(model, digits = 3, region = "NE") {
  # Obtain posterior samples
  posterior_samps <- as_draws_df(model)
  # Apply expit function to \beta_0 to get P(Y=1|X=0)
  prev <- plogis(posterior_samps$b_Intercept)
  
  # Summary of prevalence
  prevalence_summary <- data.frame(
    region = region, 
    mean = as.numeric(format(round(mean(prev), digits = digits), 
                  digits = digits, nsmall = digits)),
    lower = as.numeric(format(round(quantile(prev, 0.025), digits = digits), 
                   digits = digits, nsmall = digits)),
    upper = as.numeric(format(round(quantile(prev, 0.975), digits = digits), 
                   digits = digits, nsmall = digits)))
  
  return(list(prevalence_summary = prevalence_summary))
}


# Get region prevalence estimates for all parasites
# Use marginal models with no covariates
### (1) Ascaris
prev_reg_ascaris <- get_reg_prev(model = ascaris_bin_marg)
### (2) Hookworm
prev_reg_hook <- get_reg_prev(model = hook_bin_marg)
### (3) Trichuris
prev_reg_trich <- get_reg_prev(model = trich_bin_marg)
### (4) Strongyloides
prev_reg_strongyloides <- get_reg_prev(model = strongyloides_bin_marg)
### (5) H.nana
prev_reg_h_nana <- get_reg_prev(model = h_nana_bin_marg)
### (6) S.mansoni
prev_reg_s_mansoni <- get_reg_prev(model = s_mansoni_bin_marg)
### (7) Helms
prev_reg_helms <- get_prev_1_reg(model = helms_bin_marg)
### (8) E-coli
prev_reg_e_coli <- get_prev_1_reg(model = e_coli_bin_marg)

# Set NA for helms and e-coli regions 2-5
na_regions <- data.frame(region = c("SE", "SW", "WC", "CP"), 
                         mean = as.numeric(NA), 
                         lower = as.numeric(NA),
                         upper = as.numeric(NA))
prev_reg_helms$prevalence_summary <- rbind(prev_reg_helms$prevalence_summary, 
                                           na_regions)
prev_reg_e_coli$prevalence_summary <- rbind(prev_reg_e_coli$prevalence_summary, 
                                            na_regions)

# Append together all parasite prevalence estimates
prev_reg_all <- rbind(prev_reg_ascaris$prevalence_summary, 
                      prev_reg_hook$prevalence_summary, 
                      prev_reg_trich$prevalence_summary, 
                      prev_reg_strongyloides$prevalence_summary, 
                      prev_reg_h_nana$prevalence_summary, 
                      prev_reg_s_mansoni$prevalence_summary, 
                      prev_reg_helms$prevalence_summary, 
                      prev_reg_e_coli$prevalence_summary)
prev_reg_all$parasite <- rep(c("Ascaris", "Hookworm", "Trichuris", "Strongyloides", 
                           "H. nana", "S. mansoni", "Helms", "E. coli"), each = 5)
# Set NA for helms and e-coli regions 2-5
prev_reg_all[c(32:35, 37:40), c(2:4)] <- NA
prev_reg_all[c(32:35, 37:40), 1] <- rep(c("SE", "SW", "WC", "CP"), 2)

# Remove regions 2-5 (SE, SW, WC, CP) for helms and e-coli, since these were 
# only collected in MAHERY
# Note: for h.nana and s.mansoni, NE estimates only use data from Darwin

# # Display table
# grid::grid.newpage()
# gridExtra::grid.table(prev_reg_all_subset)
# kable(prev_reg_all_subset, format = "latex", booktabs = TRUE)

# # Save dataframe of region prevalences
# write.csv(prev_reg_all_subset, file = paste0(wd, res_dir, "Tables/prev_reg_all.csv"),
#           row.names = FALSE)


# Function to tidy up the prevalences to display in a table
tidy_prev <- function(prev_summ) {
  tidy_summ <- paste0(format(prev_summ$mean * 100, 1), " (", 
                      format(prev_summ$lower * 100, 1), ", ",
                      format(prev_summ$upper * 100, 1), ")")
  return(tidy_summ)
}

# Create nicer table
prev_reg_table <- data.frame(
  Region = c("NE", "SE", "SW", "WC", "CP"))
prev_reg_table$`A. lumbricoides` <- 
  tidy_prev(prev_reg_ascaris$prevalence_summary)
prev_reg_table$`T. trichiura` <- 
  tidy_prev(prev_reg_trich$prevalence_summary)
prev_reg_table$`Hookworm` <- 
  tidy_prev(prev_reg_hook$prevalence_summary)
prev_reg_table$`Strongyloides` <- 
  tidy_prev(prev_reg_strongyloides$prevalence_summary)
prev_reg_table$`H. nana` <- 
  tidy_prev(prev_reg_h_nana$prevalence_summary)
prev_reg_table$`S. mansoni` <- 
  tidy_prev(prev_reg_s_mansoni$prevalence_summary)
prev_reg_table$`Other Helminths` <- 
  tidy_prev(prev_reg_helms$prevalence_summary)
prev_reg_table$`E. coli` <- 
  tidy_prev(prev_reg_e_coli$prevalence_summary)

# Set NA for helms and e-coli regions 2-5
prev_reg_table[2:5, c("Other Helminths", "E. coli")] <- NA

# # Save nice table of region prevalences
# write_xlsx(prev_reg_table, 
#            path = paste0(wd, res_dir, "Tables/prev_reg_table.xlsx"))


#============== (3) Plot region-specific prevalences ===========================

# Prepare plotting data
prev_reg_all_map <- prev_reg_all %>%
  rename(Prevalence = mean,
         Region = region,
         Parasite = parasite) %>%
  mutate(Region = factor(Region, levels = c("NE", "SE", "SW", "WC", "CP")),
         Parasite = factor(Parasite, levels = c("Ascaris", "Trichuris", "Hookworm", 
                                                "Strongyloides", "H. nana", 
                                                "S. mansoni", "Helms", "E. coli"),
                           labels = c("A.lumbricoides", "T.trichiura", "Hookworm", 
                                      "Strongyloides", "H.nana", "S.mansoni", 
                                      "Helminths", "E.coli")))
parasites <- c("A.lumbricoides", "T.trichiura", "Hookworm", 
               "Strongyloides", "H.nana", "S.mansoni", 
               "Helminths", "E.coli")

### Creating heatmap of region-specific prevalences for all parasites within a 
### single facetted plot
prev_reg_all_map %>% 
  ggplot(aes(x = Parasite, y = fct_rev(Region), fill = (Prevalence * 100))) + 
  facet_grid(~ Parasite, scales = "free_x") + 
  geom_tile(color = "white",
            lwd = 1.5,
            linetype = 1) +
  scale_fill_gradient(low = "#fef0d9", high = "#d7301f", na.value = "white") +
  geom_text(aes(label = round(Prevalence * 100, 1)), 
            color = "black", size = 4) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "bottom",
        legend.direction = "horizontal",
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(face = "bold", size = 11)) + 
  guides(fill = guide_colourbar(title = "Prevalence (%)",
                                barwidth = 6,
                                barheight = 1.2)) +
  labs(y = "Region")

# # Save plot
# ggsave(filename = paste0(wd, res_dir, "Figures/region_heatmap_facetted.png"),
#        width = 11, height = 5, units = "in")


# Updated version so regions are from north to south and parasites are ordered 
# from highest to lowest prevalence
region_order <- c("NE", "CP", "SE", "WC", "SW")
parasite_order <- c("T.trichiura", "A.lumbricoides", "Helminths", "E.coli",
                    "H.nana", "Hookworm", "S.mansoni", "Strongyloides")

p_heatmap <- prev_reg_all_map %>%
  mutate(
    Region   = factor(Region, levels = region_order),
    Parasite = factor(Parasite, levels = parasite_order),
    prev_pct = Prevalence * 100
  ) %>%
  ggplot(aes(x = Parasite, y = fct_rev(Region), fill = prev_pct)) +
  facet_grid(~ Parasite, scales = "free_x") +
  geom_tile(color = "white", lwd = 1.5, linetype = 1) +
  scale_fill_gradientn(
    colours = c("#FFFEF5", "#FFF6D5", "#FEECC2", "#FED7A2", "#FDBB84", "#FC8D59", "#D7301F", "#8F477D", "#612DA4"),
    values  = scales::rescale(c(0, 0.2, 0.5, 1, 5, 20, 35, 55, 65)),
    limits  = c(0, 65),
    oob     = scales::squish,
    na.value = "white"
  ) + 
  geom_text(aes(label = round(prev_pct, 1)), color = "black", size = 4) +
  theme_bw() +
  theme(
    panel.spacing.x = unit(0.2, "lines"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "bottom",
    legend.direction = "horizontal",
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(face = "bold")
  ) +
  guides(fill = guide_colourbar(
    title = "Prevalence (%)",
    barwidth = 15,
    barheight = 1.2
  )) +
  labs(y = "Region")
p_heatmap
# ggsave(filename = paste0(wd, res_dir, "Figures/region_heatmap_facetted_new.png"),
#        width = 8.2, height = 4.8, units = "in")

# Remove helminths and e.coli from plot
p_heatmap_removed <- prev_reg_all_map %>%
  filter(Parasite != "Helminths" & Parasite != "E.coli") %>%
  mutate(
    Region   = factor(Region, levels = region_order),
    Parasite = factor(Parasite, levels = parasite_order),
    prev_pct = Prevalence * 100
  ) %>%
  ggplot(aes(x = Parasite, y = fct_rev(Region), fill = prev_pct)) +
  facet_grid(~ Parasite, scales = "free_x") +
  geom_tile(color = "white", lwd = 1.5, linetype = 1) +
  scale_fill_gradientn(
    colours = c("#FFFEF5", "#FFF6D5", "#FEECC2", "#FED7A2", "#FDBB84", "#FC8D59", "#D7301F", "#8F477D", "#612DA4"),
    values  = scales::rescale(c(0, 0.2, 0.5, 1, 5, 20, 35, 55, 65)),
    limits  = c(0, 65),
    oob     = scales::squish,
    na.value = "white"
  ) + 
  geom_text(aes(label = round(prev_pct, 1)), color = "black", size = 4) +
  theme_bw() +
  theme(
    panel.spacing.x = unit(0.5, "lines"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "bottom",
    legend.direction = "horizontal",
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(face = "bold")
  ) +
  guides(fill = guide_colourbar(
    title = "Prevalence (%)",
    barwidth = 15,
    barheight = 1.2
  )) +
  labs(y = "Region")
p_heatmap_removed
# ggsave(filename = paste0(wd, res_dir, "Figures/region_heatmap_facetted_final.png"),
#        width = 7, height = 4.8, units = "in")



### Add in map of site locations

# Read in site midpoint data
df.sites <- read_excel(
  paste0(wd, "Data/Site midpoints for ease of reference Ben 20190225c.xlsx"))

# Natural Earth Madagascar boundary
madagascar <- ne_countries(
  country = "Madagascar",
  scale = "medium",
  returnclass = "sf"
)

# Convert site coordinates to sf object
sites_sf <- st_as_sf(
  df.sites,
  coords = c("longitude", "latitude"),
  crs = 4326
)

# Manually define labeled region boxes
region_boxes <- data.frame(
  region = c("NE", "CP", "SE", "WC", "SW"),
  xmin = c(48.8, 46.2, 47.0, 43.1, 43.3),
  xmax = c(50.6, 48.0, 48.8, 44.9, 45.1),
  ymin = c(-16.2, -20.7, -22.5, -22.5, -24.3),
  ymax = c(-14.5, -19.0, -20.8, -20.8, -22.6)
)

# Label positions
region_labels <- data.frame(
  region = c("NE", "CP", "SE", "WC", "SW"),
  x = c(48.8 + 0.15, 46.2 + 0.15, 47.0 + 0.15, 44.9 - 1.1, 45.1 - 1.1),
  y = c(-14.5 - 0.15, -19.0 - 0.15, -22.5 + 0.65, -20.8 - 0.15, -24.3 + 0.65)
)

# Plot map
p.map <- ggplot() +
  geom_sf(
    data = madagascar,
    fill = "grey20",
    color = "white",
    linewidth = 0.1
  ) +
  geom_sf(
    data = sites_sf,
    color = "#72d572",
    alpha = 0.6,
    size = 2
  ) +
  geom_rect(
    data = region_boxes,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = NA,
    color = "grey80",
    linewidth = 0.6
  ) +
  geom_text(
    data = region_labels,
    aes(x = x, y = y, label = region),
    inherit.aes = FALSE,
    color = "grey80",
    fontface = "bold",
    size = 4,
    hjust = 0,
    vjust = 1
  ) +
  coord_sf(
    xlim = c(42.8, 51.2),
    ylim = c(-26, -12),
    expand = FALSE
  ) +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white", color = NA)
  )

p.map

# # Save map
# ggsave(filename = paste0(wd, res_dir, "Figures/madagascar_map.png"),
#        width = 2.5, height = 4.5, units = "in")


### Combine heatmap with madagascar map into one figure

ggarrange(p.map, p_heatmap_removed, nrow = 1, widths = c(0.25, 0.75))

# # Save combined heatmap and madagascar map
# ggsave(filename = paste0(wd, res_dir, "Figures/combined_region_heatmap_madagascar_map.png"),
#        width = 9.5, height = 4.8, units = "in")


### Creating heatmap of region-specific prevalences for all parasites, each with 
### its own legend and axes
get_reg_heatmap <- function(prev_data, parasite_index) {
  heatmap <- prev_data %>% 
    filter(Parasite == parasites[parasite_index]) %>% 
    ggplot(aes(x = Parasite, y = fct_rev(Region), fill = (Prevalence * 100))) + 
    coord_fixed() +
    geom_tile(color = "white",
              lwd = 1.5,
              linetype = 1) +
    scale_fill_gradient(low = "#fef0d9", high = "#d7301f", na.value = "white",  
                        # labels = label_number(accuracy = 2)) + 
                        labels = function(x) round(x, digits = 1)) +
    # geom_text(aes(label = format(round(Prevalence * 100, 1), 
    #                              digits = 1, nsmall = 1)), 
    #           color = "black", size = 4) +   
    geom_text(aes(label = round(Prevalence * 100, 1)), 
              color = "black", size = 4) +  
    theme_bw() +
    theme(axis.title = element_blank(), 
          # legend.position = "none"
          legend.position = "bottom", legend.direction = "horizontal") +
          # legend.text=element_text(size=8))
    guides(fill = guide_colourbar(title = "Prev\n (%)")) + 
    theme(
      axis.text.x = element_text(size = 11, face = "bold"),  # Bold parasite names
      axis.text.y = element_text(size = 11),
      legend.title = element_text(size = 7),        # Smaller legend title
      legend.text = element_text(size = 6),         # Smaller legend text
      legend.key.size = unit(0.4, "cm"),            # Smaller legend keys
      legend.key.width = unit(0.38, "cm")           # Smaller key width (for better spacing)
    ) + 
    labs(y = "Region")
  return(heatmap)
}

# shared_legend <- get_legend(
#   get_reg_heatmap(prev_data = prev_reg_all_map, parasite_index = 3))

# Get plots for all parasites and append together
p1 <- get_reg_heatmap(prev_data = prev_reg_all_map, parasite_index = 1)
p2 <- get_reg_heatmap(prev_data = prev_reg_all_map, parasite_index = 2)
p3 <- get_reg_heatmap(prev_data = prev_reg_all_map, parasite_index = 3)
p4 <- get_reg_heatmap(prev_data = prev_reg_all_map, parasite_index = 4)
p5 <- get_reg_heatmap(prev_data = prev_reg_all_map, parasite_index = 5)
p6 <- get_reg_heatmap(prev_data = prev_reg_all_map, parasite_index = 6)
p7 <- get_reg_heatmap(prev_data = prev_reg_all_map, parasite_index = 7)
p8 <- get_reg_heatmap(prev_data = prev_reg_all_map, parasite_index = 8)

p_all <- gridExtra::grid.arrange(arrangeGrob(p1, p2, p3, p4, p5, p6, p7, p8, 
                                             nrow = 1, 
                                             widths = unit(rep(1, 8), "null"),
                                             padding = unit(0.01, "cm")))

# # Save plot
# ggsave(filename = paste0(wd, res_dir, "Figures/region_heatmap_sep.png"),
#        plot = p_all, width = 11.5, height = 4, units = "in")

# Need to manually add the legend in Adobe Illustrator!



### Create plot of region prevalences and their CIs for all parasites
prev_reg_all_map %>%
  ggplot(aes(x = Prevalence*100, y = fct_rev(Region), 
             xmin = lower*100, xmax = upper*100)) + 
  geom_errorbar(col = "darkgrey") + 
  geom_point(fill = "lightblue2", pch = 21, size = 2) + 
  theme_bw() +
  facet_wrap(~ Parasite, nrow = 2, ncol = 4, scales = "fixed") + 
  theme(panel.spacing = unit(0.3, "cm")) + 
  theme(strip.background = element_rect(fill = "lightblue2")) +
  ylab("Region") + xlab("Prevalence (%)") + xlim(c(0, 100))
# # Save plot
# ggsave(filename = paste0(wd, res_dir, "Figures/prev_region_ci.png"),
#        width = 9, height = 6, units = "in")


### Look at the magnitude of variation in regional prevalence estimates
sapply(parasites, function(x) 
  max(prev_reg_all_map$Prevalence[prev_reg_all_map$Parasite == x], na.rm = TRUE) / 
    min(prev_reg_all_map$Prevalence[prev_reg_all_map$Parasite == x], na.rm = TRUE))



#============== (4) Sanity check with manual calculations ======================
# 
# # Function to calculate regional prevalence estimates for a parasite manually
# get_reg_prev_manual <- function(model, digits = 3) {
#   # Get posterior samples for intercept and region random effects
#   region_samps <- posterior_samples(model, 
#                                    pars = c("b_Intercept", "r_region"))[, 1:6]
#   
#   reg_prevs <- sapply(2:6, function(x) 
#     plogis(region_samps[, 1] + region_samps[, x]))
#   reg_prevs_summary <- apply(reg_prevs, 2, function(x) 
#     round(quantile(x, c(0.025, 0.5, 0.975)), digits = digits))
#   
#   reg_prevs_summary <- t(as.data.frame(reg_prevs_summary))
#   
#   # # Obtain region prevalence estimates mean, median, and credible interval
#   # ascaris_bin_region <- cbind(
#   #   # Mean prevalence for each region
#   #   ascaris_bin_region,
#   #   # Median and interval for each region
#   #   broom.mixed::tidyMCMC(
#   #     coda::as.mcmc(plogis(region_samps[, 1] + region_samps[, -1])),
#   #     # exp(region_samps[, 1] + region_samps[,-1]) /
#   #     #   (1 + exp(region_samps[, 1] + region_samps[, -1]))),
#   #     estimate.method = "median", conf.int = T, conf.level = c(0.9),
#   #     conf.method = 'quantile'))
#   
#   return(reg_prevs_summary)
# }
# 
# 
# # Compare prevalences with manual calculations
# prev_reg_ascaris <- get_region_prev(model = ascaris_bin_marg)
# prev_reg_ascaris_manual <- get_region_prev_manual(model = ascaris_bin_marg)
# 
# prev_reg_ascaris
# prev_reg_ascaris_manual
# 
# # Compare with crude estimates
# sapply(1:5, function(x) 
#   round(mean(fec_key_bin$ascaris_bin[fec_key_bin$region == x], na.rm = TRUE), 3))
# sapply(1:5, function(x) 
#   length(fec_key_bin$ascaris_bin[fec_key_bin$region == x]))


