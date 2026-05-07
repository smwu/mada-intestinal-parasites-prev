#===============================================================================
# Intestinal Parasites Correlations
# Author: SM Wu
# Date created: 2026/01/21
# Date updated: 2026/02/16
# Purpose: Create exploratory figures of parasite prevalence correlations
# STEPS: 
# (1) Read in data
# (2) Get prevalence correlations
# (3) Visualizations and tables
# 
# Inputs:
#   Datasets:
#       "Model_Outputs/all_fec_bin_models_tp1.RData": Model outputs for all 
#          parasites using data from the first available timepoint
#  
# Outputs (number in parentheses indicates step in which it was generated):
#   Final outputs: 
#       (3) "Model_Outputs/Figures/cooccur.png": Heatmaps of pairwise parasite co-occurrence and correlation
#       (3) "Model_Outputs/Tables/pairwise_2x2_counts_long_ordered.csv": Pairwise 2x2 tables of parasite co-occurrence


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
library(vegan)            # jaccard
library(corrplot)         # correlation
library(pheatmap)         # heatmap

### Specify directories (change to your local paths)
wd <- "~/Documents/GitHub/mada-intestinal-parasites-prev/"  # Working directory
setwd(wd)
data_dir <- "Cleaned_Data/"  # Directory with data 
code_dir <- "Code/"  # Directory with code
res_dir <- "Model_Outputs/"  # Directory to store results

# Read in dataset used
fec_key_bin <- read.csv(paste0(wd, data_dir, "fec_key_bin_tp1.csv"))


# ====== Helper functions 

# Calculate Jaccard similarity matrix
calc_jaccard <- function(df) {
  # Transpose because vegdist expects sites as rows, species as columns
  jaccard_dist <- vegdist(t(df), method = "jaccard", binary = TRUE, na.rm = TRUE)
  jaccard_sim <- 1 - as.matrix(jaccard_dist)
  return(jaccard_sim)
}

# Calculate log odds ratio matrix
calc_log_or <- function(df, label_vec = NULL) {
  n <- ncol(df)
  result <- matrix(NA, n, n)
  p_values <- matrix(NA, n, n)
  tabs <- list()
  
  colnames(result) <- rownames(result) <- colnames(df)
  colnames(p_values) <- rownames(p_values) <- colnames(df)
  
  for(i in 1:n) {
    for(j in 1:n) {
      if(i == j) {
        result[i,j] <- 0
        p_values[i,j] <- 1
      } else {
        # Use only complete cases for this pair
        complete_idx <- complete.cases(df[, c(i, j)])
        if(sum(complete_idx) > 0) {
          tab <- table(df[complete_idx, i], df[complete_idx, j])
          
          if(all(dim(tab) == c(2,2))) {
            # Add 0.5 to cells to avoid log(0) - standard continuity correction
            a <- tab[2,2] + 0.5  # both present
            b <- tab[1,2] + 0.5  # j present, i absent
            c <- tab[2,1] + 0.5  # i present, j absent
            d <- tab[1,1] + 0.5  # both absent
            
            # a <- tab[2,2] # both present
            # b <- tab[1,2]  # j present, i absent
            # c <- tab[2,1]  # i present, j absent
            # d <- tab[1,1]  # both absent
            
            result[i,j] <- log((a * d) / (b * c))
            
            # Fisher's exact test for significance
            fisher_result <- fisher.test(tab)
            p_values[i,j] <- fisher_result$p.value
          }
          
          # Add 2x2 table output
          if (i < j) {
            var_i <- colnames(df)[i]
            var_j <- colnames(df)[j]
            
            lab_i <- if (!is.null(label_vec) && var_i %in% names(label_vec)) label_vec[[var_i]] else var_i
            lab_j <- if (!is.null(label_vec) && var_j %in% names(label_vec)) label_vec[[var_j]] else var_j
            
            # set dimnames to the actual parasite names
            dimnames(tab) <- list(
              c("Absent", "Present"),
              c("Absent", "Present")
            )
            names(dimnames(tab)) <- c(lab_i, lab_j)
            
            tabs[[paste0(lab_i, " vs ", lab_j)]] <- tab
          }
        }
      }
    }
  }
  return(list(log_or = result, p_values = p_values, tables = tabs))
}


# Calculate pairwise sample sizes
calc_pairwise_n <- function(df) {
  n <- ncol(df)
  result <- matrix(NA, n, n)
  colnames(result) <- rownames(result) <- colnames(df)
  
  for(i in 1:n) {
    for(j in 1:n) {
      result[i,j] <- sum(complete.cases(df[, c(i,j)]))
    }
  }
  return(result)
}

# ======================= 2) Get prevalence correlations =======================

# Define custom labels and order
parasite_labels <- c(
  "ascaris_bin" = "A. lumbricoides",
  "trich_bin" = "T. trichiura",
  "hook_bin" = "Hookworm",
  "strongyloides_bin" = "Strongyloides",
  "h_nana_bin" = "H. nana",
  "s_mansoni_bin" = "S. mansoni",
  "helms_bin" = "Helminths",
  "e_coli_bin" = "E. coli"
)

# Desired order
parasite_order <- c("ascaris_bin", "trich_bin", "hook_bin", "strongyloides_bin",
                    "h_nana_bin", "s_mansoni_bin", "helms_bin", "e_coli_bin")

# Function to reorder and relabel matrix
reorder_and_relabel <- function(mat, order_vec, label_vec) {
  # Keep only columns that exist in the data
  order_vec <- order_vec[order_vec %in% colnames(mat)]
  
  # Reorder
  mat <- mat[order_vec, order_vec]
  
  # Relabel
  new_labels <- label_vec[colnames(mat)]
  rownames(mat) <- colnames(mat) <- new_labels
  
  return(mat)
}


# After loading your data, reorder columns first
cols_present <- parasite_order[parasite_order %in% colnames(fec_key_bin)]
fec_key_bin_var <- fec_key_bin[, cols_present]

# Calculate metrics
jaccard_sim <- calc_jaccard(fec_key_bin_var)
log_or_results <- calc_log_or(fec_key_bin_var, label_vec = parasite_labels)
sample_sizes <- calc_pairwise_n(fec_key_bin_var)

# Reorder and relabel all matrices
jaccard_sim <- reorder_and_relabel(jaccard_sim, parasite_order, parasite_labels)
log_or_results$log_or <- reorder_and_relabel(log_or_results$log_or, 
                                             parasite_order, parasite_labels)
log_or_results$p_values <- reorder_and_relabel(log_or_results$p_values, 
                                               parasite_order, parasite_labels)
sample_sizes <- reorder_and_relabel(sample_sizes, parasite_order, 
                                    parasite_labels)

# ====================== 3) Visualizations and tables ==========================

create_ggplot_heatmaps <- function(jaccard_sim, log_or_results, sample_sizes) {
  
  # Jaccard
  j_levels_y <- rownames(jaccard_sim)
  j_levels_x <- colnames(jaccard_sim)
  
  jaccard_long <- as.data.frame(jaccard_sim) %>%
    rownames_to_column("Parasite1") %>%
    pivot_longer(-Parasite1, names_to = "Parasite2", values_to = "Jaccard") %>%
    mutate(
      Parasite1 = factor(Parasite1, levels = j_levels_y),
      Parasite2 = factor(Parasite2, levels = j_levels_x),
      # diagonal -> NA (will render as gray)
      Jaccard = if_else(as.character(Parasite1) == as.character(Parasite2), NA_real_, Jaccard)
    )
  
  # LogOR + p-values
  log_levels_y <- rownames(log_or_results$log_or)
  log_levels_x <- colnames(log_or_results$log_or)
  
  log_or_long <- as.data.frame(log_or_results$log_or) %>%
    rownames_to_column("Parasite1") %>%
    pivot_longer(-Parasite1, names_to = "Parasite2", values_to = "LogOR") %>%
    mutate(
      Parasite1 = factor(Parasite1, levels = log_levels_y),
      Parasite2 = factor(Parasite2, levels = log_levels_x)
    )
  
  p_long <- as.data.frame(log_or_results$p_values) %>%
    rownames_to_column("Parasite1") %>%
    pivot_longer(-Parasite1, names_to = "Parasite2", values_to = "p_value") %>%
    mutate(
      Parasite1 = factor(Parasite1, levels = log_levels_y),
      Parasite2 = factor(Parasite2, levels = log_levels_x)
    )
  
  log_or_long <- log_or_long %>%
    left_join(p_long, by = c("Parasite1", "Parasite2")) %>%
    mutate(
      Parasite1 = factor(Parasite1, levels = log_levels_y),
      Parasite2 = factor(Parasite2, levels = log_levels_x),
      significance = case_when(
        !is.na(p_value) & p_value < 0.001 ~ "***",
        !is.na(p_value) & p_value < 0.01  ~ "**",
        !is.na(p_value) & p_value < 0.05  ~ "*",
        TRUE ~ ""
      ),
      # diagonal -> NA (will render as gray)
      LogOR = if_else(as.character(Parasite1) == as.character(Parasite2), NA_real_, LogOR),
      significance = if_else(as.character(Parasite1) == as.character(Parasite2), "", significance)
    )
  
  
  # Sample size N
  n_levels_y <- rownames(sample_sizes)
  n_levels_x <- colnames(sample_sizes)
  
  n_long <- as.data.frame(sample_sizes) %>%
    rownames_to_column("Parasite1") %>%
    pivot_longer(-Parasite1, names_to = "Parasite2", values_to = "N") %>%
    mutate(
      Parasite1 = factor(Parasite1, levels = n_levels_y),
      Parasite2 = factor(Parasite2, levels = n_levels_x),
      N = as.numeric(N),
      # diagonal -> NA (will render as gray)
      N = if_else(as.character(Parasite1) == as.character(Parasite2), NA_real_, N)
    )
  

  # Plots
  p1 <- ggplot(jaccard_long, aes(Parasite1, Parasite2, fill = Jaccard)) +
    geom_tile(color = "white") +
    geom_text(aes(label = ifelse(is.na(Jaccard), "", sprintf("%.2f", Jaccard))), 
              size = 3) +
    scale_fill_gradient(low = "white", high = "indianred3", limits = c(0, 1),
                        na.value = "grey80") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
          axis.text.y = element_text(face = "italic")) +
    labs(title = "A) Jaccard Similarity", x = "", y = "") +
    coord_fixed()
  
  p2 <- ggplot(log_or_long, aes(Parasite1, Parasite2, fill = LogOR)) +
    geom_tile(color = "white") +
    geom_text(aes(label = ifelse(is.na(LogOR), "", paste0(sprintf("%.2f", LogOR), "\n", significance))),
              size = 3) +
    scale_fill_gradient2(low = "deepskyblue4", mid = "white", high = "#F28E2B", 
                         midpoint = 0, na.value = "gray80",
                         name = "Log Odds Ratio") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
          axis.text.y = element_text(face = "italic")) +
    labs(title = "B) Log Odds Ratio", x = "", y = "") +
    coord_fixed()
  
  p3 <- ggplot(n_long, aes(Parasite1, Parasite2, fill = N)) +
    geom_tile(color = "white") +
    geom_text(aes(label = ifelse(is.na(N), "", N)), size = 3) +
    scale_fill_gradient(low = "white", high = "#796dfc",
                        na.value = "gray80") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
          axis.text.y = element_text(face = "italic")) +
    labs(title = "C) Sample Sizes", x = "", y = "") +
    coord_fixed()
  
  list(p1 = p1, p2 = p2, p3 = p3)
}


# Create heatmaps
heatmaps <- create_ggplot_heatmaps(jaccard_sim, log_or_results, sample_sizes)
# Reduce space between plots
tighten <- theme(plot.margin = margin(0, 0, 0, 0))  # top, right, bottom, left (pt)
heatmaps$p1 <- heatmaps$p1 + tighten + theme(legend.position = "bottom")
heatmaps$p2 <- heatmaps$p2 + tighten + theme(legend.position = "bottom")
heatmaps$p3 <- heatmaps$p3 + tighten + theme(legend.position = "bottom")

# Arrange all three on one page
ggarrange(
  heatmaps$p1, heatmaps$p2, heatmaps$p3,
  ncol = 3, nrow = 1, align = "hv"
)

# # Save heatmap plots
# ggsave(filename = paste0(wd, res_dir, "Figures/cooccur.png"),
#        width = 13, height = 7, units = "in")

# # Export as EPS to convert to TIFF with NAAS tool
# ggsave(filename = paste0(wd, res_dir, "Figures/Fig1.eps"),
#        width = 13, height = 7, units = "in")

# ========== Export 2x2 tables

# Counts-only long table of 2x2 data, ordered by parasite_labels
make_pairwise_2x2_long_counts_ordered <- function(df, label_vec = NULL) {
  res  <- calc_log_or(df, label_vec = label_vec)  # only using $tables
  tabs <- res$tables
  
  # Desired label order (values of the named vector)
  if (!is.null(label_vec)) {
    desired_order <- unname(label_vec)
  } else {
    # fallback: use df column names if no labels provided
    desired_order <- colnames(df)
  }
  
  out <- lapply(names(tabs), function(pair_name) {
    tab <- tabs[[pair_name]]
    
    # These are the *labels* used when you set names(dimnames(tab))
    parasite_i <- names(dimnames(tab))[1]
    parasite_j <- names(dimnames(tab))[2]
    
    # Raw 2x2 counts
    count_both_absent   <- unname(tab["Absent",  "Absent"])
    count_i_absent_j_present <- unname(tab["Absent",  "Present"])
    count_i_present_j_absent <- unname(tab["Present", "Absent"])
    count_both_present  <- unname(tab["Present", "Present"])
    n_complete_cases    <- sum(tab)
    
    data.frame(
      parasite_i = parasite_i,
      parasite_j = parasite_j,
      count_both_absent = count_both_absent,
      count_i_absent_j_present = count_i_absent_j_present,
      count_i_present_j_absent = count_i_present_j_absent,
      count_both_present = count_both_present,
      n_complete_cases = n_complete_cases,
      stringsAsFactors = FALSE
    )
  })
  
  long_df <- do.call(rbind, out)
  
  # Enforce ordering by desired_order for i, then j
  long_df$parasite_i <- factor(long_df$parasite_i, levels = desired_order)
  long_df$parasite_j <- factor(long_df$parasite_j, levels = desired_order)
  
  # Sort: first parasite_i in desired order, then parasite_j in desired order
  long_df <- long_df[order(long_df$parasite_i, long_df$parasite_j), ]
  
  # Convert back to character for nicer printing/export
  long_df$parasite_i <- as.character(long_df$parasite_i)
  long_df$parasite_j <- as.character(long_df$parasite_j)
  
  rownames(long_df) <- NULL
  long_df
}

# Run it
pairwise_counts_long <- make_pairwise_2x2_long_counts_ordered(
  fec_key_bin_var,
  label_vec = parasite_labels
)

pairwise_counts_long

# # Export 2x2 counts
# write.csv(pairwise_counts_long, 
#           file = paste0(wd, res_dir, "Tables/pairwise_2x2_counts_long_ordered.csv"), 
#                         row.names = FALSE)



#==== Summary table

create_summary_table <- function(log_or_results, sample_sizes) {
  
  n_parasites <- ncol(log_or_results$log_or)
  
  # Extract upper triangle (avoid duplicates)
  summary_df <- data.frame()
  
  for(i in 1:(n_parasites-1)) {
    for(j in (i+1):n_parasites) {
      summary_df <- rbind(summary_df, data.frame(
        Pair = paste(colnames(log_or_results$log_or)[i], 
                     colnames(log_or_results$log_or)[j], sep = " - "),
        LogOR = log_or_results$log_or[i, j],
        p_value = log_or_results$p_values[i, j],
        N = sample_sizes[i, j]
      ))
    }
  }
  
  # Add interpretation
  summary_df <- summary_df %>%
    mutate(
      Association = case_when(
        p_value >= 0.05 ~ "Not significant",
        LogOR > 0 ~ "Positive co-occurrence",
        LogOR < 0 ~ "Negative co-occurrence"
      ),
      Significance = case_when(
        p_value < 0.001 ~ "***",
        p_value < 0.01 ~ "**",
        p_value < 0.05 ~ "*",
        TRUE ~ "ns"
      )
    ) %>%
    arrange(p_value)
  
  return(summary_df)
}

summary_table <- create_summary_table(log_or_results, sample_sizes)
print(summary_table)
