#===============================================================================
# Intestinal Parasites EDA and Prevalence Estimates
# Author: SM Wu
# Date created: 2025/01/28
# Date updated: 2025/11/17
# Purpose: Create exploratory figures and parasite prevalence estimates 
# STEPS: 
# (1) Read in data 
# (2) Subset to first observed timepoint
# (3) Exploratory histograms and descriptive tables
# (4) Prevalence of binary parasite presence
# (5) Display overall prevalence results
# 
# Inputs:
#   Datasets:
#       "Cleaned_Data/mah_dar_crs_fec_cleaned_avg_wide_drop_miss.csv": Cleaned 
#       parasite data in wide format, excluding those missing age or sex
#       "Cleaned_Data/mah_dar_crs_fec_cleaned_avg.csv" Cleaned fecal sample 
#       averaged data w/ age and sex for MAH, DAR, and CRS, in long format and 
#       with missingness
#  
# Outputs (number in parentheses indicates step in which it was generated):
#   Final outputs: 
#       (3) "Model_Outputs/Figures/hist.png": Histogram distributions of key parasites
#       (3) "Model_Outputs/Tables/Table1.doc": Table of sample characteristics
#       (4) "Model_Outputs/all_fec_bin_models_tp1.Rdata": All prevalence models for all parasites, using first observed timepoint
#       (5) "Model_Outputs/prev_all.csv": Dataframe of overall prevalences for all parasites
#       (5) "Model_Outputs/Figures/prev_overall_ci.png": Plot of overall prevalences and 95% credible intervals for all parasites

#============== (1) Read in data ===============================================
# Clear memory
rm(list = ls())

# Load libraries
library(tidyverse)        # for tidy data routines
library(magrittr)         # for efficient piping
library(ggplot2)          # for plotting
library(ggbeeswarm)       # for beeswarm plots
library(brms)             # for Bayesian hierarchical modeling
library(knitr)            # for displaying tables
library(tidybayes)        # for Bayesian modeling helper functions
library(grid)             # for displaying tables
library(gridExtra)        # for displaying tables
library(gtsummary)        # create descriptive table
library(ggpubr)           # displaying multiple plots
library(naniar)           # check missingness
library(loo)              # comparing models


## use all cores for speed
num_cores <- parallel::detectCores()
options(mc.cores = num_cores)

### Specify directories (change to your local paths)
# setwd("/Users/Stephanie/Documents/GitHub/Intestinal_Parasites")
wd <- "/Users/Stephanie/Documents/GitHub/Intestinal_Parasites/"  # Working directory
data_dir <- "Cleaned_Data/"  # Directory with data 
code_dir <- "Code/"  # Directory with code
res_dir <- "Model_Outputs/"  # Directory to store results

# Read in cleaned averaged parasite data in wide format
fec_all_wide_clean <- read.csv(
  paste0(wd, data_dir, "mah_dar_crs_fec_cleaned_avg_wide_drop_miss.csv"))
fec_all_wide_clean <- fec_all_wide_clean %>%
  rename(h_nana = `h..nana`,
         s_mansoni = `s..mansoni`)

# Summary of values for each parasite
summary(fec_all_wide_clean %>% 
          select(ascaris, trich, hook, strongyloides, h_nana, s_mansoni, 
                 helms, e_coli, fungus, insect, pollen))

# Read in cleaned averaged parasite data including missing age and sex
fec_all_clean_miss <- read.csv(
  paste0(wd, data_dir, "mah_dar_crs_fec_cleaned_avg.csv"))

#============== (2) Subset to first observed timepoint =========================

# Total 4443 observations: 3689 tp1, 618 tp2, 136 tp3
nrow(fec_all_wide_clean)
table(fec_all_wide_clean$time_point)
# How many observations occur at tp 2 or 3: 618 tp2, 136 tp3
mult_tp <- fec_all_wide_clean %>%
  filter(time_point > 1)
table(mult_tp$time_point, useNA = "always")

# Convert to long format
fec_all_clean <- fec_all_wide_clean %>%
  pivot_longer(cols = ascaris:s_mansoni, 
               names_to = "parasite", values_to = "value")

# Differences between timepoints
fec_all_clean_indiv <- fec_all_clean %>%
  group_by(project_name, region, village, household, individual, parasite) %>%
  summarise(num_tp = length(value), 
            values = list(value),
            max = max(value), min = min(value), diff = max - min, 
            mean = mean(value), sd = sd(value)) %>%
  filter(num_tp > 1)
fec_all_clean_indiv_wide <- fec_all_clean_indiv %>%
  select(-c(max, min, mean, sd)) %>%
  pivot_wider(names_from = parasite, values_from = c(diff, values))
# Number of individuals w/ more than one tp: 569 (only 1 with 3 tps)
table(fec_all_clean_indiv_wide$num_tp)
summary(fec_all_clean_indiv$sd)

## For each parasite, subset to first timepoint with data for each individual
# This is done because individual-level random intercepts are unstable (not 
# enough repeated data)
fec_all <- fec_all_clean %>%
  arrange(project_name, region, village, household, individual, parasite, time_point) %>%
  group_by(project_name, region, village, household, individual, parasite, ) %>%
  mutate(first_tp = first(na.omit(time_point))) %>%
  filter(time_point == first_tp) %>% ungroup()
# Wide format
fec_all_wide <- fec_all %>%
  pivot_wider(names_from = parasite, values_from = value)
# Number of individuals: 3872 (was 3901 before remove age/sex)
nrow(fec_all_wide)
# Individuals by dataset: mahery 597, darwin 760, crs 2515
table(fec_all_wide$project_name)
# Individuals by region
table(fec_all_wide$region)
# Number of households: 1035
nrow(fec_all_wide %>% select(project_name, region, village, household) %>% 
       distinct())



# # Save cross-sectional intestinal parasite data subsetted to first available
# # timepoint for each individual
# write.csv(fec_all_wide, 
#           file = paste0(wd, data_dir, "fec_all_cross_sec.csv"), 
#           row.names = FALSE)



# Apply the same subsetting to first available timepoint for those with missing 
# age and sex
fec_all_wide_miss <- fec_all_clean_miss %>%
  arrange(project_name, region, village, household, individual, parasite, time_point) %>%
  group_by(project_name, region, village, household, individual, parasite, ) %>%
  mutate(first_tp = first(na.omit(time_point))) %>%
  filter(time_point == first_tp) %>% ungroup() %>%
  pivot_wider(names_from = parasite, values_from = value)  %>%
  rename(h_nana = `h. nana`,
         s_mansoni = `s. mansoni`)
# Number of individuals: 3902
nrow(fec_all_wide_miss)
# Individuals by dataset
table(fec_all_wide_miss$project_name)
# Individuals by region
table(fec_all_wide_miss$region)
# Number of households: 1041
nrow(fec_all_wide_miss %>% select(project_name, region, village, household) %>% 
       distinct())


#======== (3) Exploratory histograms and descriptive tables ====================



fec_all %>% 
  ggplot(aes(x = parasite, y = value, fill = parasite, col = parasite)) + 
  geom_violin(width = 1, trim = TRUE, alpha = 0.1,
              linewidth = 0.1) +
  geom_quasirandom(dodge.width = 0.7, cex=0.15, alpha=0.5, method = "tukeyDense",
                   width = 0.25) +
  #geom_boxplot(width = 0.8, position = dodge) + 
  stat_summary(fun=mean, geom="point", size=1, col = "black") + 
  theme_bw() + 
  facet_wrap(~parasite, scales = "free", ncol = 3) + 
  theme(legend.position = "right",
        strip.background = element_rect(fill="aliceblue")) + 
  coord_flip()

# Histogram of values for each parasite
fec_all %>% 
  mutate(parasite = factor(parasite, 
                          levels = c("ascaris", "trich", "hook", "strongyloides", 
                                     "h_nana", "s_mansoni", "helms", "e_coli", 
                                     "fungus", "insect", "pollen"),
                          labels = c("A.lumbricoides", "T.trichiura", "Hookworm", 
                                     "Strongyloides", "H.nana", "S.mansoni", 
                                     "Helminths", "E.coli", "Fungus", "Insect", 
                                     "Pollen"))) %>%
  ggplot(aes(x = value, fill = parasite)) + 
  theme_bw() + 
  geom_histogram() + 
  facet_wrap(~parasite, scales = "free", ncol = 4) + 
  theme(legend.position = "none")

# Histogram for ones where may want to use a different cutoff
fec_all %>% 
  filter(parasite %in% c("ascaris", "trich", "hook", "strongyloides", 
                         "h_nana", "s_mansoni", "helms", "e_coli")) %>%
  mutate(parasite = factor(parasite, 
                           levels = c("ascaris", "trich", "hook", "strongyloides", 
                                      "h_nana", "s_mansoni", "helms", "e_coli"),
                           labels = c("A.lumbricoides", "T.trichiura", "Hookworm", 
                                      "Strongyloides", "H.nana", "S.mansoni", 
                                      "Helminths", "E.coli"))) %>%
  ggplot(aes(x = value, fill = parasite)) + 
  theme_bw() + 
  geom_histogram() + 
  facet_wrap(~parasite, scales = "free", ncol = 4) + 
  theme(legend.position = "none",
        strip.text = element_text(size = 14)) + 
  xlab("Egg Count") + ylab("Frequency")
# # Save plot
# ggsave(filename = paste0(res_dir, "Figures/hist.png"), width = 12, height = 6,
#        units = "in")


# Ascaris with custom bins
binned_data <- cut(unlist(fec_all %>% filter(parasite == "ascaris") %>% 
                            select(value)), 
                   breaks = c(0, 0.25, 1, 5, 10, 50, 100, 500, 1100), 
                   include.lowest = TRUE, right = FALSE)
table(binned_data)


## Create descriptive table

# Rename columns
table_data <- fec_all_wide %>%
  ungroup() %>%
  rename(`A.lumbricoides` = ascaris,
         `T.trichiura` = trich,
         Hookworm = hook,
         Strongyloides = strongyloides,
         `H.nana` = `h_nana`,
         `S.mansoni` = `s_mansoni`,
         `Helminths` = helms,
         `E.coli` = e_coli)
# table_data <- fec_all_wide %>%
#   ungroup() %>%
#   rename(`A.lumbricoides` = ascaris,
#          `T.trichiura` = trich,
#          Hookworm = hook,
#          Strongyloides = strongyloides,
#          `H.nana` = `h_nana`,
#          `S.mansoni` = `s_mansoni`,
#          `Helminths` = helms,
#          `E.coli` = e_coli)

parasite_vars <- c("A.lumbricoides", "T.trichiura", "Hookworm", "Strongyloides", 
                   "H.nana", "S.mansoni", "Helminths", "E.coli")

# Create binary variables (any/none), but preserve NA for missing
df <- table_data %>%
  mutate(across(
    all_of(parasite_vars),
    ~ case_when(
      is.na(.) ~ NA_integer_,
      . > 0 ~ 1L,
      TRUE ~ 0L
    ),
    .names = "{.col}_bin"
  ))

# Make binary parasite variables into factors with "Yes", "No", and NA
df <- df %>%
  mutate(across(
    ends_with("_bin"),
    ~ factor(., levels = c(1, 0), labels = c("Yes", "No"))
  ))

# Combine all needed variables
all_vars <- c("ageyears", "age_cat", "sex", "A.lumbricoides", "A.lumbricoides_bin", 
              "T.trichiura", "T.trichiura_bin", "Hookworm", "Hookworm_bin",
              "Strongyloides", "Strongyloides_bin", "H.nana", "H.nana_bin", 
              "S.mansoni", "S.mansoni_bin", "Helminths", "Helminths_bin", 
              "E.coli", "E.coli_bin")

# Define which are continuous and which are categorical
continuous_vars <- parasite_vars  # All original parasite counts
categorical_vars <- c("age_cat", "sex", paste0(parasite_vars, "_bin"))

# Ensure region is a factor
df <- df %>% mutate(region = as.factor(region),
                    across(all_of(parasite_vars), as.numeric))

# Create Table 1
table1 <- df %>%
  select(all_of(all_vars), region) %>%
  tbl_summary(
    by = region,
    missing = "ifany",  # Display missingness if it exists
    type = list(
      ageyears ~ "continuous",
      all_of(continuous_vars) ~ "continuous",
      # Helminths ~ "continuous",
      # Strongyloides ~ "continuous",
      all_of(categorical_vars) ~ "categorical"
    ),
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = list(
      # Global defaults
      all_continuous() ~ c(4, 4),  # 4 decimals for mean, 4 for SD
      all_categorical() ~ c(0, 1), # 0 decimals for counts, 1 for percentages
      
      # Variable-specific digits
      ageyears ~ c(1, 1),        # mean 1 decimal, SD 1 decimal
      A.lumbricoides ~ c(1, 1),  # mean 1 decimal, SD 1 decimal
      T.trichiura ~ c(2, 2),     # 2 decimals mean, 2 SD
      Hookworm ~ c(2, 2),        # 2 decimals mean, 2 SD
      Helminths ~ c(2, 2),       # 2 decimals mean, 2 SD
      E.coli ~ c(2, 2)           # 2 decimals mean, 2 SD
    ),
        
    label = list(
      ageyears ~ "Age (years)",
      age_cat ~ "Age group",
      sex ~ "Sex",
      `A.lumbricoides` ~ "A. lumbricoides (mean)",
      `A.lumbricoides_bin` ~ "A. lumbricoides (any)",
      `T.trichiura` ~ "T. trichiura (mean)",
      `T.trichiura_bin` ~ "T. trichiura (any)",
      Hookworm ~ "Hookworm (mean)",
      Hookworm_bin ~ "Hookworm (any)",
      Strongyloides ~ "Strongyloides (mean)",
      Strongyloides_bin ~ "Strongyloides (any)",
      `H.nana` ~ "H. nana (mean)",
      `H.nana_bin` ~ "H. nana (any)",
      `S.mansoni` ~ "S. mansoni (mean)",
      `S.mansoni_bin` ~ "S. mansoni (any)",
      Helminths ~ "Helminths (mean)",
      Helminths_bin ~ "Helminths (any)",
      `E.coli` ~ "E. coli (mean)",
      `E.coli_bin` ~ "E. coli (any)"
    )
  ) %>%
  add_overall(last = TRUE) %>%
  modify_header(label ~ "**Variable**")

# Export
export_table <- table1 %>%
  as_flex_table() %>%
  flextable::width(width = 1)
flextable::save_as_docx(
  export_table, 
  path = paste0(wd, res_dir, "Tables/Table1.docx"), width = 1)



# Check missingness
miss_df <- fec_all_wide_miss %>% select(-c(pollen, insect, fungus))%>%
  rename(`A.lumbricoides` = ascaris,
         `T.trichiura` = trich,
         Hookworm = hook,
         Strongyloides = strongyloides,
         `H.nana` = `h_nana`,
         `S.mansoni` = `s_mansoni`,
         `Helminths` = helms,
         `E.coli` = e_coli) %>%
  mutate(across(
    all_of(parasite_vars),
    ~ case_when(
      is.na(.) ~ NA_integer_,
      . > 0 ~ 1L,
      TRUE ~ 0L
    ),
    .names = "{.col}_bin"
  ))
gg_miss_var(miss_df)
vis_miss(miss_df)

missing_summary <- miss_df %>%
  summarise(across(everything(), ~sum(is.na(.)))) %>%
  pivot_longer(everything(), names_to = "Variable", values_to = "Missing_N") %>%
  mutate(Missing_Perc = round(100 * Missing_N / nrow(df), 1))
missing_summary

# Create Table 1
table1_miss <- miss_df %>%
  select(all_of(all_vars), region) %>%
  tbl_summary(
    by = region,
    missing = "ifany",  # Display missingness if it exists
    type = list(
      ageyears ~ "continuous",
      all_of(continuous_vars) ~ "continuous",
      # Helminths ~ "continuous",
      # Strongyloides ~ "continuous",
      all_of(categorical_vars) ~ "categorical"
    ),
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = list(
      all_continuous() ~ c(1, 1),  # 1 decimal for mean and SD
      all_categorical() ~ c(0, 1)  # 0 decimals for counts, 1 for percentages
    ),
    
    label = list(
      ageyears ~ "Age (years)",
      age_cat ~ "Age group",
      sex ~ "Sex",
      `A.lumbricoides` ~ "A. lumbricoides (mean)",
      `A.lumbricoides_bin` ~ "A. lumbricoides (any)",
      `T.trichiura` ~ "T. trichiura (mean)",
      `T.trichiura_bin` ~ "T. trichiura (any)",
      Hookworm ~ "Hookworm (mean)",
      Hookworm_bin ~ "Hookworm (any)",
      Strongyloides ~ "Strongyloides (mean)",
      Strongyloides_bin ~ "Strongyloides (any)",
      `H.nana` ~ "H. nana (mean)",
      `H.nana_bin` ~ "H. nana (any)",
      `S.mansoni` ~ "S. mansoni (mean)",
      `S.mansoni_bin` ~ "S. mansoni (any)",
      Helminths ~ "Helminths (mean)",
      Helminths_bin ~ "Helminths (any)",
      `E.coli` ~ "E. coli (mean)",
      `E.coli_bin` ~ "E. coli (any)"
    )
  ) %>%
  add_overall(last = TRUE) %>%
  modify_header(label ~ "**Variable**")
table1_miss

#========= (4) Prevalence of binary parasite presence ==========================

# Focus on key parasites (drop insects, fungus, pollen)
fec_key <- fec_all_wide %>%
  select(-c(fungus, insect, pollen))

# Create binary variables
fec_key_bin <- fec_key %>%
  mutate(across(ascaris:trich, ~ as.integer(. > 0), 
                .names = "{.col}_bin"))

# Set default for sex and age categories
fec_key_bin <- fec_key_bin %>%
  # Reference category is set to largest categories: females, age [20,50) 
  mutate(sex = factor(sex, levels = c(1, 0), labels = c("F", "M")),
         age_cat = factor(age_cat, levels = c("[20,50)", "[0,2)", "[2,5)",  
                                              "[5,12)", "[12,20)", "[50,Inf)")))

# # Save dataset with binary variables used for model fitting
# write.csv(fec_key_bin, file = paste0(wd, data_dir, "fec_key_bin_tp1.csv"),
#           row.names = FALSE)

# Function to obtain prevalence estimates from brms model
get_prev <- function(model, digits = 3) {
  # Obtain posterior samples
  posterior_samps <- as_draws_df(model)
  # Apply expit function to \beta_0 to get P(Y=1|X=0)
  prev <- plogis(posterior_samps$b_Intercept)
  # Summary of prevalence
  prev_summary <- c(mean = mean(prev), sd = sd(prev), 
                    quantile(prev, probs = c(0.025, 0.6, 0.975)))
  prev_summary <- format(round(prev_summary, digits), digits = digits,
                         nsmall = digits)
  # Convert prevalences to numeric
  prev_summary <- as.numeric(prev_summary)
  
  return(prev_summary)
}


#============== Run Bayesian models for prevalence estimation ==================

## (1) Ascaris
set.seed(1)
# model with just main effects
ascaris_bin_main <- brm(ascaris_bin ~ sex + age_cat + 
                          (1|region/village/household), 
                        data = fec_key_bin, 
                        family = bernoulli(), 
                        control = list(adapt_delta = 0.99, max_treedepth = 15),
                        iter = 4000, warmup = 1000, thin = 3,
                        chains = 4, cores = num_cores)
set.seed(1)
# model with no covariates
ascaris_bin_marg <- brm(ascaris_bin ~ (1 | region/village/household), 
                        data = fec_key_bin, 
                        family = bernoulli(), 
                        control = list(adapt_delta = 0.99),
                        iter = 4000, warmup = 1000, thin = 3,
                        chains = 4, cores = num_cores)
set.seed(1)
# model with age-sex interactions
ascaris_bin_int <- brm(ascaris_bin ~ sex * age_cat + 
                          (1|region/village/household), 
                        data = fec_key_bin, 
                        family = bernoulli(), 
                        control = list(adapt_delta = 0.99),
                        iter = 4000, warmup = 1000, thin = 3,
                        chains = 4, cores = num_cores)
# # Test with no covariates or clustering
# ascaris_bin_test <- brm(ascaris_bin ~ 1, 
#                         data = fec_key_bin, 
#                         family = bernoulli(), 
#                         control = list(adapt_delta = 0.99), 
#                         iter=3000, warmup = 800, thin = 3,
#                         chains = 4, cores = num_cores)
# get_prev(ascaris_bin_test)
# # Check crude prevalence
# mean(fec_key_bin$ascaris_bin, na.rm = TRUE)
# sd(fec_key_bin$ascaris_bin, na.rm = TRUE)


## (2) Hookworm
set.seed(2)
# model with just main effects
hook_bin_main <- brm(hook_bin ~ sex + age_cat + 
                       (1|region/village/household), 
                     data = fec_key_bin, 
                     family = bernoulli(), 
                     control = list(adapt_delta = 0.99),
                     iter = 4000, warmup = 1000, thin = 3,
                     chains = 4, cores = num_cores) #1div
set.seed(2)
# model with no covariates
hook_bin_marg <- brm(hook_bin ~ (1 | region/village/household), 
                     data = fec_key_bin, 
                     family = bernoulli(), 
                     control = list(adapt_delta = 0.999, max_treedepth = 15),
                     iter = 4000, warmup = 1000, thin = 3,
                     chains = 4, cores = num_cores) 
set.seed(2)
# model with age-sex interactions
hook_bin_int <- brm(hook_bin ~ sex * age_cat + 
                       (1|region/village/household), 
                     data = fec_key_bin, 
                     family = bernoulli(), 
                     control = list(adapt_delta = 0.99),
                     iter = 4000, warmup = 1000, thin = 3,
                     chains = 4, cores = num_cores) #1div


## (3) Trichuris
set.seed(3)
# model with just main effects
trich_bin_main <- brm(trich_bin ~ sex + age_cat + 
                        (1|region/village/household), 
                      data = fec_key_bin, 
                      family = bernoulli(), 
                      control = list(adapt_delta = 0.99, max_treedepth = 15),
                      iter = 4000, warmup = 1000, thin = 3,
                      chains = 4, cores = num_cores)  
set.seed(3)
# model with no covariates
trich_bin_marg <- brm(trich_bin ~ (1 | region/village/household), 
                      data = fec_key_bin, 
                      family = bernoulli(), 
                      control = list(adapt_delta = 0.999, max_treedepth = 15),
                      iter = 4000, warmup = 1000, thin = 3,
                      chains = 4, cores = num_cores) 
set.seed(3)
# model with age-sex interactions
trich_bin_int <- brm(trich_bin ~ sex * age_cat + 
                        (1|region/village/household), 
                      data = fec_key_bin, 
                      family = bernoulli(), 
                      control = list(adapt_delta = 0.99, max_treedepth = 15),
                      iter = 4000, warmup = 1000, thin = 3,
                      chains = 4, cores = num_cores) 


## (4) Strongyloides: n = 2266
set.seed(4)
# model with just main effects
strongyloides_bin_main <- brm(strongyloides_bin ~ sex + age_cat + 
                                (1|region/village/household), 
                              data = fec_key_bin, 
                              family = bernoulli(), 
                              control = list(adapt_delta = 0.995),
                              iter = 4000, warmup = 1000, thin = 3,
                              chains = 4, cores = num_cores) #1div
set.seed(4)
# model with no covariates
strongyloides_bin_marg <- brm(strongyloides_bin ~ 
                                (1 | region/village/household), 
                              data = fec_key_bin, 
                              family = bernoulli(), 
                              control = list(adapt_delta = 0.99),
                              iter = 4000, warmup = 1000, thin = 3,
                              chains = 4, cores = num_cores) 
set.seed(4)
# model with age-sex interactions: 
strongyloides_bin_int <- brm(strongyloides_bin ~ sex * age_cat + 
                                (1|region/village/household), 
                              data = fec_key_bin, 
                              family = bernoulli(), 
                              control = list(adapt_delta = 0.999, 
                                             max_treedepth = 15),
                              iter = 4000, warmup = 1000, thin = 3,
                              chains = 4, cores = num_cores) 


## (5) H.nana (not measured in MAHERY)
set.seed(5)
# model with just main effects
h_nana_bin_main <- brm(h_nana_bin ~ sex + age_cat + 
                         (1|region/village/household), 
                       data = fec_key_bin, 
                       family = bernoulli(), 
                       control = list(adapt_delta = 0.99),
                       iter = 4000, warmup = 1000, thin = 3,
                       chains = 4, cores = num_cores) 
set.seed(5)
# model with no covariates
h_nana_bin_marg <- brm(h_nana_bin ~ (1 | region/village/household), 
                       data = fec_key_bin, 
                       family = bernoulli(), 
                       control = list(adapt_delta = 0.99),
                       iter = 4000, warmup = 1000, thin = 3,
                       chains = 4, cores = num_cores) #1div
set.seed(5)
# model with age-sex interactions 
h_nana_bin_int <- brm(h_nana_bin ~ sex * age_cat + 
                         (1|region/village/household), 
                       data = fec_key_bin, 
                       family = bernoulli(), 
                       control = list(adapt_delta = 0.99, max_treedepth = 15),
                       iter = 4000, warmup = 1000, thin = 3,
                       chains = 4, cores = num_cores) 


## (6) S.mansoni (not measured in MAHERY)
set.seed(6)
# model with just main effects
s_mansoni_bin_main <- brm(s_mansoni_bin ~ sex + age_cat + 
                            (1|region/village/household), 
                          data = fec_key_bin, 
                          family = bernoulli(), 
                          control = list(adapt_delta = 0.99),
                          iter = 4000, warmup = 1000, thin = 3,
                          chains = 4, cores = num_cores) 
set.seed(6)
# model with no covariates
s_mansoni_bin_marg <- brm(s_mansoni_bin ~ (1 | region/village/household), 
                          data = fec_key_bin, 
                          family = bernoulli(), 
                          control = list(adapt_delta = 0.99),
                          iter = 4000, warmup = 1000, thin = 3,
                          chains = 4, cores = num_cores) 
set.seed(6)
# model with age-sex interactions 
s_mansoni_bin_int <- brm(s_mansoni_bin ~ sex * age_cat + 
                            (1|region/village/household), 
                          data = fec_key_bin, 
                          family = bernoulli(), 
                          control = list(adapt_delta = 0.99, max_treedepth = 15),
                          iter = 4000, warmup = 1000, thin = 3,
                          chains = 4, cores = num_cores) 


## (7) Helminths (not measured in Darwin or CRS)
set.seed(7)
# model with just main effects
helms_bin_main <- brm(helms_bin ~ sex + age_cat + 
                        (1|region/village/household), 
                      data = fec_key_bin, 
                      family = bernoulli(), 
                      control = list(adapt_delta = 0.99, max_treedepth = 15),
                      iter = 4000, warmup = 1000, thin = 3,
                      chains = 4, cores = num_cores) #1div
set.seed(7)
# model with no covariates
helms_bin_marg <- brm(helms_bin ~ (1 | region/village/household), 
                      data = fec_key_bin, 
                      family = bernoulli(), 
                      control = list(adapt_delta = 0.995, max_treedepth = 15),
                      iter = 4000, warmup = 1000, thin = 3,
                      chains = 4, cores = num_cores) 
set.seed(7)
# model with age-sex interactions
helms_bin_int <- brm(helms_bin ~ sex * age_cat + 
                        (1|region/village/household), 
                      data = fec_key_bin, 
                      family = bernoulli(), 
                      control = list(adapt_delta = 0.999, max_treedepth = 15),
                      iter = 4000, warmup = 1000, thin = 3,
                      chains = 4, cores = num_cores) 


## (8) E-coli (not measured in Darwin or CRS)
set.seed(8)
# model with just main effects
e_coli_bin_main <- brm(e_coli_bin ~ sex + age_cat + 
                         (1|region/village/household), 
                       data = fec_key_bin, 
                       family = bernoulli(), 
                       control = list(adapt_delta = 0.999, max_treedepth = 15),
                       iter = 4000, warmup = 1000, thin = 3,
                       chains = 4, cores = num_cores)   
set.seed(8)
# model with no covariates
e_coli_bin_marg <- brm(e_coli_bin ~ (1 | region/village/household), 
                       data = fec_key_bin, 
                       family = bernoulli(), 
                       control = list(adapt_delta = 0.9995, max_treedepth = 15),
                       iter = 4000, warmup = 1000, thin = 3,
                       chains = 4, cores = num_cores) #1div 
set.seed(8)
# model with age-sex interactions
e_coli_bin_int <- brm(e_coli_bin ~ sex * age_cat + 
                         (1|region/village/household), 
                       data = fec_key_bin, 
                       family = bernoulli(), 
                       control = list(adapt_delta = 0.999, max_treedepth = 15),
                       iter = 4000, warmup = 1000, thin = 3,
                       chains = 4, cores = num_cores) 


# Examine prevalences
get_prev(ascaris_bin_main) # n=3872
get_prev(ascaris_bin_int) # n=3872
get_prev(ascaris_bin_marg) # n=3901 (old version. Now all n=3872)
summary(fec_key_bin$ascaris_bin)
get_prev(hook_bin_main) # n=3872
get_prev(hook_bin_int) # n=3872
get_prev(hook_bin_marg) # n=3872
summary(fec_key_bin$hook_bin)
get_prev(trich_bin_main) # n=3872
get_prev(trich_bin_int) # n=3872
get_prev(trich_bin_marg) # n=3872
summary(fec_key_bin$trich_bin)
get_prev(strongyloides_bin_main) # n=3872
get_prev(strongyloides_bin_int) # n=3872
get_prev(strongyloides_bin_marg) # n=3872
summary(fec_key_bin$strongyloides_bin)
get_prev(h_nana_bin_main) # n=3275
get_prev(h_nana_bin_int) # n=3275
get_prev(h_nana_bin_marg) # n=3275
summary(fec_key_bin$h_nana_bin)
get_prev(s_mansoni_bin_main) # n=3275
get_prev(s_mansoni_bin_int) # n=3275
get_prev(s_mansoni_bin_marg) # n=3275
summary(fec_key_bin$s_mansoni_bin)
get_prev(helms_bin_main) # n=597
get_prev(helms_bin_int) # n=597
get_prev(helms_bin_marg) # n=603 (old version. Now all 597)
summary(fec_key_bin$helms_bin)
get_prev(e_coli_bin_main) # n=597
get_prev(e_coli_bin_int) # n=597
get_prev(e_coli_bin_marg) # n=97
summary(fec_key_bin$e_coli_bin)

# Compare models with leave-one-out cross validation
loo_compare(loo(ascaris_bin_int), loo(ascaris_bin_main)) # main better
loo_compare(loo(trich_bin_int), loo(trich_bin_main)) # main better
loo_compare(loo(hook_bin_int), loo(hook_bin_main)) # main better
loo_compare(loo(strongyloides_bin_int), loo(strongyloides_bin_main)) # main better
loo_compare(loo(h_nana_bin_int), loo(h_nana_bin_main)) # interaction better
loo_compare(loo(s_mansoni_bin_int), loo(s_mansoni_bin_main)) # main better
loo_compare(loo(helms_bin_int), loo(helms_bin_main)) # main better
loo_compare(loo(e_coli_bin_int), loo(e_coli_bin_main)) # main better



# # Save all models using data from the first available timepoint
save(ascaris_bin_main, ascaris_bin_marg, ascaris_bin_int,
     trich_bin_main, trich_bin_marg, trich_bin_int,
     hook_bin_main, hook_bin_marg, hook_bin_int,
     strongyloides_bin_main, strongyloides_bin_marg, strongyloides_bin_int,
     h_nana_bin_main, h_nana_bin_marg, h_nana_bin_int,
     s_mansoni_bin_main, s_mansoni_bin_marg, s_mansoni_bin_int,
     helms_bin_main, helms_bin_marg, helms_bin_int,
     e_coli_bin_main, e_coli_bin_marg, e_coli_bin_int,
     file = paste0(wd, res_dir, "all_fec_bin_models_tp1.Rdata"))

# # Load models
# load(paste0(wd, res_dir, "all_fec_bin_models_tp1.Rdata"))

# Model diagnostics for interaction models
g1 <- pp_check(ascaris_bin_int)
g2 <- pp_check(trich_bin_int)
g3 <- pp_check(hook_bin_int)
g4 <- pp_check(strongyloides_bin_int)
g5 <- pp_check(h_nana_bin_int)
g6 <- pp_check(s_mansoni_bin_int)
g7 <- pp_check(helms_bin_int)
g8 <- pp_check(e_coli_bin_int)
ggarrange(plotlist = list(g1, g2, g3, g4, g5, g6, g7, g8), nrow = 3, ncol = 3, 
          labels = c("A. lumbricoides", "T. trichiura", "Hookworm", 
                      "Strongyloides", "H. nana", "S. mansoni", 
                       "Helminths", "E. coli"),
          font.label = list(size = 11, color = "black"))

# # Save plot
# ggsave(filename = paste0(wd, res_dir, "Figures/pp_check_int.png"),
#        width = 7, height = 5, units = "in")


# Model diagnostics for marginal models
g1 <- pp_check(ascaris_bin_marg)
g2 <- pp_check(trich_bin_marg)
g3 <- pp_check(hook_bin_marg)
g4 <- pp_check(strongyloides_bin_marg)
g5 <- pp_check(h_nana_bin_marg)
g6 <- pp_check(s_mansoni_bin_marg)
g7 <- pp_check(helms_bin_marg)
g8 <- pp_check(e_coli_bin_marg)
ggarrange(plotlist = list(g1, g2, g3, g4, g5, g6, g7, g8), nrow = 3, ncol = 3, 
          labels = c("A. lumbricoides", "T. trichiura", "Hookworm", 
                     "Strongyloides", "H. nana", "S. mansoni", 
                     "Helminths", "E. coli"),
          # common.legend = TRUE,
          font.label = list(size = 11, color = "black"))

# # Save plot
# ggsave(filename = paste0(wd, res_dir, "Figures/pp_check_marg.png"),
#        width = 7, height = 5, units = "in")

#============= (5) Display overall prevalence results ==========================

# Initialize data frame of results 
prev_df <- as.data.frame(matrix(NA, nrow = 8, ncol = 6))
prev_df[, 1] <- c("A. lumbricoides", "T. trichiura", "Hookworm", "Strongyloides", 
                  "H. nana", "S. mansoni", "Helminths", "E. coli")
colnames(prev_df) <- c("Parasite", "Mean", "SD", "2.5%", "Median", "97.5%")
prev_df[1, -1] <- get_prev(model = ascaris_bin_marg) * 100
prev_df[2, -1] <- get_prev(model = trich_bin_marg) * 100
prev_df[3, -1] <- get_prev(model = hook_bin_marg) * 100
prev_df[4, -1] <- get_prev(model = strongyloides_bin_marg) * 100
prev_df[5, -1] <- get_prev(model = h_nana_bin_marg) * 100
prev_df[6, -1] <- get_prev(model = s_mansoni_bin_marg) * 100
prev_df[7, -1] <- get_prev(model = helms_bin_marg) * 100
prev_df[8, -1] <- get_prev(model = e_coli_bin_marg) * 100
# Reorganize columns
prev_df <- prev_df[, c("Parasite", "Mean", "Median", "2.5%", "97.5%")]

# Display table
kable(prev_df, format = "latex", booktabs = TRUE)
grid::grid.newpage()
gridExtra::grid.table(prev_df)

# # Save dataframe of prevalences
# write.csv(prev_df, file = paste0(wd, res_dir, "Tables/prev_all.csv"),
#           row.names = FALSE)


### Create plot of region prevalences and their CIs for all parasites
prev_df %>%
  rename(Prevalence = Mean,
         lower = `2.5%`,
         upper = `97.5%`) %>%
  ggplot(aes(x = Prevalence, y = Parasite, xmin = lower, xmax = upper)) + 
  geom_errorbar(col = "darkgrey") + 
  geom_point(aes(x = Prevalence, shape = "Mean"), size = 2, fill = "cyan3") + 
  geom_point(aes(x = Median, shape = "Median"), color = "orange3", size = 2) +
  scale_shape_manual(values = c("Mean" = 21, "Median" = 4)) +  # Custom shapes
  labs(shape = "Summary") + 
  theme_bw() +
  ylab("Parasite") + xlab("Prevalence (%)") 
# # Save plot
# ggsave(filename = paste0(wd, res_dir, "Figures/prev_overall_ci.png"),
#        width = 5, height = 4.5, units = "in")

# ### Create halfeye plot of prevalence distributions
# # Get posterior samples of prevalences for all parasites
# all_samps <- as.data.frame(
#   sapply(list(ascaris_bin_marg, hook_bin_marg, trich_bin_marg, 
#               strongyloides_bin_marg, h_nana_bin_marg, s_mansoni_bin_marg, 
#               helms_bin_marg, e_coli_bin_marg), 
#          function(x) 
#   plogis(as_draws_df(x, variable = "b_Intercept")$b_Intercept)))
# colnames(all_samps) <- c("Ascaris", "Hookworm", "Trichuris", "Strongyloides", 
#                          "H. nana", "S. mansoni", "Helms", "E. coli")
# all_samps %>% 
#   pivot_longer(cols = everything(), names_to = "Parasite", 
#                values_to = "Prevalence") %>%
#   ggplot(aes(x = Prevalence)) +
#   facet_wrap(~ Parasite, scales = "free_x") + 
#   stat_halfeye() +
#   labs(x = "Prevalence", y = NULL) +
#   theme_bw() +
#   theme(legend.position = "right")


#============== Sanity checks ==================================================

# ## Sanity check of overall prevalence obtained from main model
# # Get the empirical distribution of sex and age_cat
# weighting_table <- fec_key_bin %>%
#   drop_na(ascaris_bin, sex, age_cat) %>%
#   group_by(sex, age_cat) %>%
#   summarise(weight = n(), .groups = "drop") %>%
#   mutate(weight = weight / sum(weight))  # Convert counts to proportions
# # Create new data with unique (sex, age_cat) combinations
# new_data <- weighting_table %>%
#   select(sex, age_cat) %>%
#   mutate( # ensure factor levels match
#     sex = factor(sex, levels = levels(fec_key_bin$sex)),
#     age_cat = factor(age_cat, levels = levels(fec_key_bin$age_cat))
#   )
# # Generate posterior draws of the linear predictor
# posterior_samples <- posterior_epred(ascaris_bin_main, 
#                                      newdata = new_data, 
#                                      allow_new_levels = TRUE)
# # Compute the weighted prevalence for each posterior draw
# weighted_prevalence <- apply(posterior_samples, 1, function(p) {
#   sum(p * weighting_table$weight)  # Weighted sum of probabilities
# })
# # Summarize results
# overall_prevalence <- mean(weighted_prevalence)
# credible_interval <- quantile(weighted_prevalence, probs = c(0.025, 0.975))
# # Print results
# cat("Overall Prevalence Estimate:", overall_prevalence, "\n")
# cat("95% Credible Interval:", credible_interval, "\n")
# get_prev(ascaris_bin_main) # n=2266
# get_prev(ascaris_bin_marg) # n=3708
# summary(fec_key_bin$ascaris_bin)

#============== [OLD] Prevalence models with individual random intercepts ======

# # Re-define data to include all timepoints
# 
# # Focus on key parasites (drop insects, fungus, pollen)
# fec_key_all <- fec_all_wide_clean %>%
#   select(-c(fungus, insect, pollen))
# 
# # Create binary variables
# fec_key_bin <- fec_key_clean %>%
#   rename(h_nana = h_nana, 
#          s_mansoni = s_mansoni) %>%
#   mutate(across(ascaris:`s_mansoni`, ~ as.integer(. > 0), 
#                 .names = "{.col}_bin"))
# 
# # Set default for sex and age categories
# fec_key_bin <- fec_key_bin %>%
#   # Reference category is set to largest categories: females, age [20,50) 
#   mutate(sex = factor(sex, levels = c(1, 0), labels = c("F", "M")),
#          age_cat = factor(age_cat, levels = c("[20,50)", "[0,2)", "[2,5)",  
#                                               "[5,12)", "[12,20)", "[50,Inf)")))
# 
# 
# ## (1) Ascaris
# set.seed(1)
# # model with just main effects
# ascaris_bin_main <- brm(ascaris_bin ~ sex + age_cat + 
#                           (1|region/village/household/individual), 
#                         data = fec_key_bin, 
#                         family = bernoulli(), 
#                         # control = list(adapt_delta = 0.99), 
#                         iter = 4000, warmup = 1000, thin = 3,
#                         chains = 4, cores = num_cores)
# set.seed(1)
# # model with no covariates
# ascaris_bin_marg <- brm(ascaris_bin ~ (1 | region/village/household/individual), 
#                         data = fec_key_bin, 
#                         family = bernoulli(), 
#                         # control = list(adapt_delta = 0.99), 
#                         iter = 4000, warmup = 1000, thin = 3,
#                         chains = 4, cores = num_cores)
# # # Test with no covariates or clustering
# # ascaris_bin_test <- brm(ascaris_bin ~ 1, 
# #                         data = fec_key_bin, 
# #                         family = bernoulli(), 
# #                         control = list(adapt_delta = 0.99), 
# #                         iter=3000, warmup = 800, thin = 3,
# #                         chains = 4, cores = num_cores)
# # get_prev(ascaris_bin_test)
# # # Check crude prevalence
# # mean(fec_key_bin$ascaris_bin, na.rm = TRUE)
# # sd(fec_key_bin$ascaris_bin, na.rm = TRUE)
# 
# 
# ## (2) Hookworm
# set.seed(2)
# # model with just main effects
# hook_bin_main <- brm(hook_bin ~ sex + age_cat + 
#                        (1|region/village/household/individual), 
#                      data = fec_key_bin, 
#                      family = bernoulli(), 
#                      # control = list(adapt_delta = 0.99), 
#                      iter = 4000, warmup = 1000, thin = 3,
#                      chains = 4, cores = num_cores)
# 
# # model with no covariates
# hook_bin_marg <- brm(hook_bin ~ (1 | region/village/household/individual), 
#                      data = fec_key_bin, 
#                      family = bernoulli(), 
#                      # control = list(adapt_delta = 0.99), 
#                      iter = 4000, warmup = 1000, thin = 3,
#                      chains = 4, cores = num_cores)
# 
# 
# ## (3) Trichuris
# set.seed(3)
# # model with just main effects
# trich_bin_main <- brm(trich_bin ~ sex + age_cat + 
#                         (1|region/village/household/individual), 
#                       data = fec_key_bin, 
#                       family = bernoulli(), 
#                       # control = list(adapt_delta = 0.99), 
#                       iter = 4000, warmup = 1000, thin = 3,
#                       chains = 4, cores = num_cores)
# 
# # model with no covariates
# trich_bin_marg <- brm(trich_bin ~ (1 | region/village/household/individual), 
#                       data = fec_key_bin, 
#                       family = bernoulli(), 
#                       # control = list(adapt_delta = 0.99), 
#                       iter = 4000, warmup = 2000, thin = 3,
#                       chains = 4, cores = num_cores)
# 
# 
# ## (4) Strongyloides
# set.seed(4)
# # model with just main effects
# strongyloides_bin_main <- brm(strongyloides_bin ~ sex + age_cat + 
#                        (1|region/village/household/individual), 
#                      data = fec_key_bin, 
#                      family = bernoulli(), 
#                      # control = list(adapt_delta = 0.99), 
#                      iter = 4000, warmup = 1000, thin = 3,
#                      chains = 4, cores = num_cores)
# 
# # model with no covariates
# strongyloides_bin_marg <- brm(strongyloides_bin ~ 
#                                 (1 | region/village/household/individual), 
#                      data = fec_key_bin, 
#                      family = bernoulli(), 
#                      # control = list(adapt_delta = 0.99), 
#                      iter = 4000, warmup = 1000, thin = 3,
#                      chains = 4, cores = num_cores)
# 
# 
# ## (5) H.nana
# set.seed(5)
# # model with just main effects
# h_nana_bin_main <- brm(h_nana_bin ~ sex + age_cat + 
#                        (1|region/village/household/individual), 
#                      data = fec_key_bin, 
#                      family = bernoulli(), 
#                      # control = list(adapt_delta = 0.99), 
#                      iter = 4000, warmup = 1000, thin = 3,
#                      chains = 4, cores = num_cores)
# 
# # model with no covariates
# h_nana_bin_marg <- brm(h_nana_bin ~ (1 | region/village/household/individual), 
#                      data = fec_key_bin, 
#                      family = bernoulli(), 
#                      # control = list(adapt_delta = 0.99), 
#                      iter = 4000, warmup = 1000, thin = 3,
#                      chains = 4, cores = num_cores)
# 
# 
# ## (6) S.mansoni
# set.seed(6)
# # model with just main effects
# s_mansoni_bin_main <- brm(s_mansoni_bin ~ sex + age_cat + 
#                         (1|region/village/household/individual), 
#                       data = fec_key_bin, 
#                       family = bernoulli(), 
#                       # control = list(adapt_delta = 0.99), 
#                       iter = 4000, warmup = 1000, thin = 3,
#                       chains = 4, cores = num_cores)
# 
# # model with no covariates
# s_mansoni_bin_marg <- brm(s_mansoni_bin ~ (1 | region/village/household/individual), 
#                       data = fec_key_bin, 
#                       family = bernoulli(), 
#                       # control = list(adapt_delta = 0.99), 
#                       iter = 4000, warmup = 1000, thin = 3,
#                       chains = 4, cores = num_cores)
# 
# 
# ## (7) Helminths
# set.seed(7)
# # model with just main effects
# helms_bin_main <- brm(helms_bin ~ sex + age_cat + 
#                        (1|region/village/household/individual), 
#                      data = fec_key_bin, 
#                      family = bernoulli(), 
#                      # control = list(adapt_delta = 0.99), 
#                      iter = 4000, warmup = 1000, thin = 3,
#                      chains = 4, cores = num_cores)
# 
# # model with no covariates
# helms_bin_marg <- brm(helms_bin ~ (1 | region/village/household/individual), 
#                      data = fec_key_bin, 
#                      family = bernoulli(), 
#                      # control = list(adapt_delta = 0.99), 
#                      iter = 4000, warmup = 1000, thin = 3,
#                      chains = 4, cores = num_cores)
# 
# 
# ## (8) E-coli
# set.seed(8)
# # model with just main effects
# e_coli_bin_main <- brm(e_coli_bin ~ sex + age_cat + 
#                        (1|region/village/household/individual), 
#                      data = fec_key_bin, 
#                      family = bernoulli(), 
#                      # control = list(adapt_delta = 0.99), 
#                      iter = 4000, warmup = 1000, thin = 3,
#                      chains = 4, cores = num_cores)
# 
# # model with no covariates
# e_coli_bin_marg <- brm(e_coli_bin ~ (1 | region/village/household/individual), 
#                      data = fec_key_bin, 
#                      family = bernoulli(), 
#                      # control = list(adapt_delta = 0.99), 
#                      iter = 4000, warmup = 1000, thin = 3,
#                      chains = 4, cores = num_cores)
# 
# 
# # Examine prevalences
# get_prev(ascaris_bin_main)
# get_prev(ascaris_bin_marg) # check
# summary(fec_key_bin$ascaris_bin)
# get_prev(trich_bin_main) # check
# get_prev(trich_bin_marg) # check
# summary(fec_key_bin$trich_bin)
# get_prev(hook_bin_main)
# get_prev(hook_bin_marg) 
# summary(fec_key_bin$hook_bin)
# get_prev(strongyloides_bin_main)
# get_prev(strongyloides_bin_marg)
# summary(fec_key_bin$strongyloides_bin)
# get_prev(h_nana_bin_main)
# get_prev(h_nana_bin_marg) # check
# summary(fec_key_bin$h_nana_bin)
# get_prev(s_mansoni_bin_main)
# get_prev(s_mansoni_bin_marg) # check
# summary(fec_key_bin$s_mansoni_bin)
# get_prev(helms_bin_main) # check
# get_prev(helms_bin_marg) # check
# summary(fec_key_bin$helms_bin)
# get_prev(e_coli_bin_main) # check
# get_prev(e_coli_bin_marg) # check
# summary(fec_key_bin$e_coli_bin)


# # Save all models
# save(ascaris_bin_main, ascaris_bin_marg,
#      trich_bin_main, trich_bin_marg,
#      hook_bin_main, hook_bin_marg,
#      strongyloides_bin_main, strongyloides_bin_marg,
#      h_nana_bin_main, h_nana_bin_marg,
#      s_mansoni_bin_main, s_mansoni_bin_marg,
#      helms_bin_main, helms_bin_marg,
#      e_coli_bin_main, e_coli_bin_marg,
#      file = paste0(wd, res_dir, "all_fec_bin_models.Rdata"))

