################################################################################
# Admixture / distance plots for snps  Identified by Outflank, Bayescan, 
# and PCAadapt
# 
# This script loads SNPs identified as outliers from three different selection 
# detection methods (Outflank, Bayescan, and PCAadapt), runs LEA for med ancestry
# and plots ancestry vs geographic distance from med.

# Author: Martin Taylor (martin.taylor@uea.ac.uk)
################################################################################

#load libraries
library(LEA).      # Model based clustering
library(dplyr)     # Data manipulation
library(dartR)     # SNP genotype handling
library(plotrix)   # Plotting utilitieslibrary(tidyr)
library(ggplot2)   # Data visualization
library(reshape2)  # Matrix manipulation
library(tidyr)     # Data transformation

# ---------------------------------------------------------------
# Load Data
# ---------------------------------------------------------------

# Load outliers detected by Bayescan
bayescan_outliers <- gl.load("output/bayescan/outlying_loci_bayescan.gl")

# Load outliers detected by Outflank
outflank_outliers <- gl.load("output/outflank/sp_outlying_loci_outflank.gl")

# Load outliers detected by PCAadapt
pcaadapt_outliers <- gl.load("output/pca_adapt/sp_outlying_loci_pca_adapt.gl")

#Load outliers detected by all 3 methods
all_outliers <- gl.load("output/outflank/sp_outlying_loci_all.gl")

# Load pairwise sea distances (precomputed)
dist_mat_sp <- as.matrix(read.table("output/sea_distances/Contour_spawning_incl_med_18_pops.txt", sep = ","))
dist.df_sp <- reshape2::melt(dist_mat_sp)  # Convert matrix to data frame

# Extract distances from Catalonia & clean column names
filtered_dist <- dist.df_sp %>% filter(Var1 == "Catalonia") %>%
  rename(origin = Var1, ices_rect = Var2, dist = value) %>%
  mutate(ices_rect = gsub("X", "", as.character(ices_rect)))

# Load ICES area data
ices_tab_uniq<-read.csv("data/ices_area_rectangle.csv", sep=",")

# check pops are correct across datasets
# Extract population levels from each dataset
pop_levels_all <- levels(all_outliers$pop)
pop_levels_pca <- levels(pcaadapt_outliers$pop)
pop_levels_out <- levels(outflank_outliers$pop)
pop_levels_bayes <- levels(bayescan_outliers$pop)

# Check if all are identical
same_levels <- identical(pop_levels_all, pop_levels_pca) &&
  identical(pop_levels_all, pop_levels_out) &&
  identical(pop_levels_all, pop_levels_bayes)

# Print result
if (same_levels) {
  print("All datasets have the same population levels.")
} else {
  print("Population levels differ across datasets.")
  
  # Optional: Show differences if they exist
  print("Differences between datasets:")
  print(list(
    "All vs PCAadapt" = all.equal(pop_levels_all, pop_levels_pca),
    "All vs Outflank" = all.equal(pop_levels_all, pop_levels_out),
    "All vs Bayescan" = all.equal(pop_levels_all, pop_levels_bayes)
  ))
}

# ---------------------------------------------------------------
# Convert Genlight Objects to GENO Files for LEA
# ---------------------------------------------------------------

#make function to convert gl to geno for LEA
save_outliers_as_geno <- function(genlight_object, output_file) {
  
  # Convert genlight object to data frame
  outliers_df <- as.data.frame(genlight_object)
  
  # Replace NA values with 9
  outliers_df[is.na(outliers_df)] <- 9
  
  # Write the transformed data to a .geno file
  write.geno(outliers_df, output_file)
  
  # Print confirmation message
  cat("File successfully saved to:", output_file, "\n")
}

#save outliers as geno files
#all outliers
save_outliers_as_geno(all_outliers, "data/LEA_input/outlier_snps_all_df.geno")
#bayescan outliers
save_outliers_as_geno(bayescan_outliers, "data/LEA_input/bayescan_snps_df.geno")
#outflank outliers
save_outliers_as_geno(outflank_outliers, "data/LEA_input/outflank_snps_df.geno")
#pcaadapt outliers
save_outliers_as_geno(pcaadapt_outliers, "data/LEA_input/pcaadapt_snps_df.geno")

# ---------------------------------------------------------------
# Run Admixture Analysis Using snmf
# ---------------------------------------------------------------

#set up snmf project runs for each file

outlier_snps_all_lea <- snmf("data/LEA_input/outlier_snps_all_df.geno",
                          K = 1:3,
                          entropy = TRUE,
                          repetitions = 10,
                          project = "new",
                          CPU = 4)

outlier_snps_bayes_lea <- snmf("data/LEA_input/bayescan_snps_df.geno",
                             K = 1:3,
                             entropy = TRUE,
                             repetitions = 10,
                             project = "new",
                             CPU = 4)

outlier_snps_outflank_lea <- snmf("data/LEA_input/outflank_snps_df.geno",
                               K = 1:3,
                               entropy = TRUE,
                               repetitions = 10,
                               project = "new",
                               CPU = 4)

outlier_snps_pcaadapt_lea <- snmf("data/LEA_input/pcaadapt_snps_df.geno",
                                  K = 1:3,
                                  entropy = TRUE,
                                  repetitions = 10,
                                  project = "new",
                                  CPU = 4)

# ---------------------------------------------------------------
# Extract and Process Q-Matrices
# ---------------------------------------------------------------

#function to extract info
extract_best_q_matrix <- function(lea_object, K = 2) {
  # Load necessary library
  library(LEA)
  
  # Identify the best run using cross-entropy
  best_run <- which.min(cross.entropy(lea_object, K = K))
  
  # Extract the Q matrix for the best run
  q_matrix <- LEA::Q(lea_object, K = K, run = best_run)
  
  # Rename columns to reflect population proportions
  colnames(q_matrix) <- paste0("P", 1:K)
  
  return(q_matrix)
}


#run function on lea outputs for each outlier method
q_mat_all <- extract_best_q_matrix(outlier_snps_all_lea, K = 2)
q_mat_bayes <- extract_best_q_matrix(outlier_snps_bayes_lea, K = 2)
q_mat_outflank <- extract_best_q_matrix(outlier_snps_outflank_lea, K = 2)
q_mat_pcaadapt <- extract_best_q_matrix(outlier_snps_pcaadapt_lea, K = 2)

#convert to df
#for all outliers
q_df_all <- q_mat_all %>% 
  as_tibble() %>% 
  # add pop and indiv info for plotting
  mutate(individual=outliers_all.gl$ind.names,
         ices_rect=outliers_all.gl$pop
  )

#for bayescan
q_df_bayes <- q_mat_bayes %>% 
  as_tibble() %>% 
  # add pop and indiv info for plotting
  mutate(individual=bayescan_outliers$ind.names,
         ices_rect=bayescan_outliers$pop
  )

#for outflank
q_df_outflank <- q_mat_outflank %>% 
  as_tibble() %>% 
  # add pop and indiv info for plotting
  mutate(individual=outflank_outliers$ind.names,
         ices_rect=outflank_outliers$pop
  )

#for pcaadapt
q_df_pcaadapt <- q_mat_pcaadapt %>% 
  as_tibble() %>% 
  # add pop and indiv info for plotting
  mutate(individual=pcaadapt_outliers$ind.names,
         ices_rect=pcaadapt_outliers$pop
  )

# ---------------------------------------------------------------
# Summarize Admixture Proportions (Mean, SD, SE)
# ---------------------------------------------------------------

summary_all <- q_df_all %>%
  group_by(ices_rect) %>%
  summarise(across(where(is.numeric), 
                   list(mean = mean, sd = sd, se = ~ sd(.) / sqrt(n())), 
                   .names = "{.col}_{.fn}"))

summary_bayes <- q_df_bayes %>%
  group_by(ices_rect) %>%
  summarise(across(where(is.numeric), 
                   list(mean = mean, sd = sd, se = ~ sd(.) / sqrt(n())), 
                   .names = "{.col}_{.fn}"))

summary_outflank <- q_df_outflank %>%
  group_by(ices_rect) %>%
  summarise(across(where(is.numeric), 
                   list(mean = mean, sd = sd, se = ~ sd(.) / sqrt(n())), 
                   .names = "{.col}_{.fn}"))

summary_pcaadapt <- q_df_pcaadapt %>%
  group_by(ices_rect) %>%
  summarise(across(where(is.numeric), 
                   list(mean = mean, sd = sd, se = ~ sd(.) / sqrt(n())), 
                   .names = "{.col}_{.fn}"))


# Merge summary data with ICES areas and distances
summary_join_all <- summary_all %>%
  left_join(ices_tab_uniq, by = "ices_rect") %>%
  left_join(filtered_dist, by = "ices_rect") %>%
  mutate(dist2 = dist - 1859,  # Adjust distance from Portugal_N
         group = ifelse(ices_rect %in% c("Catalonia", "Portugal_N", "Portugal_S"), 1, 0))

summary_join_bayes <- summary_bayes %>%
  left_join(ices_tab_uniq, by = "ices_rect") %>%
  left_join(filtered_dist, by = "ices_rect") %>%
  mutate(dist2 = dist - 1859,  # Adjust distance from Portugal_N
         group = ifelse(ices_rect %in% c("Catalonia", "Portugal_N", "Portugal_S"), 1, 0))

summary_join_pcaadapt <- summary_pcaadapt %>%
  left_join(ices_tab_uniq, by = "ices_rect") %>%
  left_join(filtered_dist, by = "ices_rect") %>%
  mutate(dist2 = dist - 1859,  # Adjust distance from Portugal_N
         group = ifelse(ices_rect %in% c("Catalonia", "Portugal_N", "Portugal_S"), 1, 0))

summary_join_outflank <- summary_outflank %>%
  left_join(ices_tab_uniq, by = "ices_rect") %>%
  left_join(filtered_dist, by = "ices_rect") %>%
  mutate(dist2 = dist - 1859,  # Adjust distance from Portugal_N
         group = ifelse(ices_rect %in% c("Catalonia", "Portugal_N", "Portugal_S"), 1, 0))

# ---------------------------------------------------------------
# Plot Admixture Proportions Against Distance
# ---------------------------------------------------------------

#all outliers
all_zoom <- ggplot(summary_join_all[c(1:15), ], aes(x = dist2, y = P2_mean, color = ices_area)) +
  geom_errorbar(aes(ymin = P2_mean - (2 * P2_se), ymax = P2_mean + (2 * P2_se)), width = 0.1) +
  geom_smooth(aes(group = 1), method = "lm", se = FALSE, colour = "black", size = 0.5) +
  geom_point(size = 3) +
  scale_color_manual(values = c("4.c" = "#1D72F5", "7.a" = "#DF0101", "7.b" = "#77CE61", 
                                "7.d" = "#FF9326", "7.e" = "#A945FF", "7.f" = "#0089B2", 
                                "MED" = "#FDF060", "Portugal_N" = "#FFA6B2", "Portugal_S" = "#BFF217")) +
  theme_bw() +
  labs(x = "Distance from Portugal_N (km)", y = "Admixture proportion") +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15, margin = margin(t = 20)))

#bayes outliers
bayes_zoom <- ggplot(summary_join_bayes[c(1:15), ], aes(x = dist2, y = P2_mean, color = ices_area)) +
  geom_errorbar(aes(ymin = P2_mean - (2 * P2_se), ymax = P2_mean + (2 * P2_se)), width = 0.1) +
  geom_smooth(aes(group = 1), method = "lm", se = FALSE, colour = "black", size = 0.5) +
  geom_point(size = 3) +
  scale_color_manual(values = c("4.c" = "#1D72F5", "7.a" = "#DF0101", "7.b" = "#77CE61", 
                                "7.d" = "#FF9326", "7.e" = "#A945FF", "7.f" = "#0089B2", 
                                "MED" = "#FDF060", "Portugal_N" = "#FFA6B2", "Portugal_S" = "#BFF217")) +
  theme_bw() +
  labs(x = "Distance from Portugal_N (km)", y = "Admixture proportion") +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15, margin = margin(t = 20)))

#outflank outliers - p1 and p2 reversed here
outflank_zoom <- ggplot(summary_join_outflank[c(1:15), ], aes(x = dist2, y = P1_mean, color = ices_area)) +
  geom_errorbar(aes(ymin = P1_mean - (2 * P1_se), ymax = P1_mean + (2 * P1_se)), width = 0.1) +
  geom_smooth(aes(group = 1), method = "lm", se = FALSE, colour = "black", size = 0.5) +
  geom_point(size = 3) +
  scale_color_manual(values = c("4.c" = "#1D72F5", "7.a" = "#DF0101", "7.b" = "#77CE61", 
                                "7.d" = "#FF9326", "7.e" = "#A945FF", "7.f" = "#0089B2", 
                                "MED" = "#FDF060", "Portugal_N" = "#FFA6B2", "Portugal_S" = "#BFF217")) +
  theme_bw() +
  labs(x = "Distance from Portugal_N (km)", y = "Admixture proportion") +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15, margin = margin(t = 20)))

#pcaadapt outliers
pcaadapt_zoom <- ggplot(summary_join_pcaadapt[c(1:15), ], aes(x = dist2, y = P2_mean, color = ices_area)) +
  geom_errorbar(aes(ymin = P2_mean - (2 * P2_se), ymax = P2_mean + (2 * P2_se)), width = 0.1) +
  geom_smooth(aes(group = 1), method = "lm", se = FALSE, colour = "black", size = 0.5) +
  geom_point(size = 3) +
  scale_color_manual(values = c("4.c" = "#1D72F5", "7.a" = "#DF0101", "7.b" = "#77CE61", 
                                "7.d" = "#FF9326", "7.e" = "#A945FF", "7.f" = "#0089B2", 
                                "MED" = "#FDF060", "Portugal_N" = "#FFA6B2", "Portugal_S" = "#BFF217")) +
  theme_bw() +
  labs(x = "Distance from Portugal_N (km)", y = "Admixture proportion") +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15, margin = margin(t = 20)))


#save plot
pdf("output/admixture/outliers_dist_all.pdf")
all_zoom
dev.off()

pdf("output/admixture/outliers_dist_bayes.pdf")
bayes_zoom
dev.off()

pdf("output/admixture/outliers_dist_outflank.pdf")
outflank_zoom
dev.off()

pdf("output/admixture/outliers_dist_pcaadapt.pdf")
pcaadapt_zoom
dev.off()

# Fit linear model & print equation

#all outliers
fit.all <- lm(P2_mean ~ dist2, data = summary_join_all[c(1:15), ])
summary(fit.all)

#bayes
fit.bayes <- lm(P2_mean ~ dist2, data = summary_join_bayes[c(1:15), ])
summary(fit.bayes)

#outflank
fit.outflank <- lm(P2_mean ~ dist2, data = summary_join_outflank[c(1:15), ])
summary(fit.outflank)

#pcaadapt
fit.pcaadapt <- lm(P2_mean ~ dist2, data = summary_join_pcaadapt[c(1:15), ])
summary(fit.pcaadapt)

# Function to extract model statistics
extract_lm_stats <- function(data, formula, subset_rows = 1:15) {
  # Fit linear model
  fit <- lm(formula, data = data[subset_rows, ])
  
  # Extract summary statistics
  fit_summary <- summary(fit)
  r_squared <- fit_summary$r.squared
  p_value <- coef(fit_summary)[2, 4]  # Extract p-value for slope
  
  # Extract coefficients and format equation
  coefficients <- fit$coefficients
  equation <- paste0("Y = ", round(coefficients[1], 2), 
                     " + ", round(coefficients[2], 2), " * dist2 + e")
  
  # Return results as a named vector
  return(c(R2 = round(r_squared, 4), 
           p_value = round(p_value, 4), 
           equation = equation))
}

# Create a data frame with results for all models
model_stats <- data.frame(
  Model = c("All Outliers", "Bayes", "Outflank", "PCAadapt"),
  do.call(rbind, list(
    extract_lm_stats(summary_join_all, P2_mean ~ dist2),
    extract_lm_stats(summary_join_bayes, P2_mean ~ dist2),
    extract_lm_stats(summary_join_outflank, P1_mean ~ dist2),#other way round for outflank
    extract_lm_stats(summary_join_pcaadapt, P2_mean ~ dist2)
  ))
)

# Print the table
print(model_stats)

#Model           R2      p_value      equation
#1 All Outliers 0.0193  0.6213 Y = 0 + 0 * dist2 + e
#2        Bayes 0.0504  0.4211 Y = 0 + 0 * dist2 + e
#3     Outflank  7e-04  0.9261 Y = 0 + 0 * dist2 + e
#4     PCAadapt 0.0435  0.4556 Y = 0 + 0 * dist2 + e
