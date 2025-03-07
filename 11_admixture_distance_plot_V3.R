
######################################################
#Script to output admixture vs distance from MED plots 
#Fig 7 sea bass manuscript
#
##admixture proportion with geographic distance
#med ancestry with distance from Catalonia plot
#
#
#Martin Taylor (martin.taylor@uea.ac.uk)
######################################################


# Load libraries
library(LEA)       # Model-based clustering
library(dplyr)     # Data manipulation
library(dartR)     # Genetic analysis
library(plotrix)   # Additional plotting functions
library(ggplot2)   # Data visualization
library(reshape2)  # Reshape matrices

# Load pairwise sea distances (precomputed)
dist_mat_sp <- as.matrix(read.table("output/sea_distances/Contour_spawning_incl_med_18_pops.txt", sep = ","))
dist.df_sp <- reshape2::melt(dist_mat_sp)  # Convert matrix to data frame

# Extract distances from Catalonia & clean column names
filtered_dist <- dist.df_sp %>% filter(Var1 == "Catalonia") %>%
  rename(origin = Var1, ices_rect = Var2, dist = value) %>%
  mutate(ices_rect = gsub("X", "", as.character(ices_rect)))

# Load spawning + Mediterranean dataset (41022 loci, 343 individuals)
bass_reduced_ices_area_spawning_MED <- gl.load("data/spawning_incl_med.gl")

# Convert to data frame & replace NAs with 9
bass_reduced_ices_area_spawning_MED_df <- as.data.frame(bass_reduced_ices_area_spawning_MED)
bass_reduced_ices_area_spawning_MED_df[is.na(bass_reduced_ices_area_spawning_MED_df)] <- 9

# Save genotype file for snmf analysis
write.geno(bass_reduced_ices_area_spawning_MED_df, "data/LEA_input/bass_reduced_ices_area_spawning_MED_df.geno")

# Run snmf clustering (K = 1 to 3)
spawning_incl_med_K2_3 <- snmf("data/LEA_input/bass_reduced_ices_area_spawning_MED_df.geno",
                               K = 1:3, entropy = TRUE, repetitions = 10, 
                               project = "new", CPU = 4)

# Plot cross-entropy to identify the best K
plot(spawning_incl_med_K2_3, col = "blue", pch = 19, cex = 1.2)

# Select best K = 2 run
best_sp_med <- which.min(cross.entropy(spawning_incl_med_K2_3, K = 2))

# Extract Q matrix & name columns
q_mat_sp_med <- as.data.frame(LEA::Q(spawning_incl_med_K2_3, K = 2, run = best_sp_med))
colnames(q_mat_sp_med) <- paste0("P", 1:2)
q_mat_sp_med$ices_rect <- bass_reduced_ices_area_spawning_MED$pop

# Summarize admixture proportions (mean, SD, SE)
summary <- q_mat_sp_med %>%
  group_by(ices_rect) %>%
  summarise(across(where(is.numeric), 
                   list(mean = mean, sd = sd, se = ~ sd(.) / sqrt(n())), 
                   .names = "{.col}_{.fn}"))

# Load ICES area data
ices_tab_uniq<-read.csv("data/ices_area_rectangle.csv", sep=",")

# Merge summary data with ICES areas and distances
summary_join <- summary %>%
  left_join(ices_tab_uniq, by = "ices_rect") %>%
  left_join(filtered_dist, by = "ices_rect") %>%
  mutate(dist2 = dist - 1859,  # Adjust distance from Portugal_N
         group = ifelse(ices_rect %in% c("Catalonia", "Portugal_N", "Portugal_S"), 1, 0))

# Plot admixture proportions against distance
spawn_zoom <- ggplot(summary_join[c(1:15), ], aes(x = dist2, y = P2_mean, color = ices_area)) +
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
pdf("output/admixture/spawning_dist_zoom_sept24.pdf")
spawn_zoom
dev.off()


# Fit linear model & print equation
fit.lm <- lm(P2_mean ~ dist2, data = summary_join[c(1:15), ])
summary(fit.lm)

cc <- fit.lm$coefficients
eqn <- paste0("Y = ", round(cc[1], 2), " + ", round(cc[2], 2), " * dist2 + e")
print(eqn)
