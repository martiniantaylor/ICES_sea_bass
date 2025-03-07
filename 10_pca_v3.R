# PCA Analysis for Figure 5

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(adegenet)
library(dartR)

# Load spawning data for PCA (including Atlantic, Mediterranean, and Portugal)
pca_data <- gl.load("data/spawning_reduced_incl_portN_portS_med")

# Perform PCA with 20 factors
pca_med <- glPca(pca_data, nf = 20)

# Compute variance explained
eig_total_med <- sum(pca_med$eig)
variance_explained <- function(index) {
  formatC((pca_med$eig[index] / eig_total_med) * 100, format = "f", digits = 2)
}

PC1_variance <- variance_explained(1)
PC2_variance <- variance_explained(2)

# Convert PCA scores to dataframe and add population labels
pca_scores <- as.data.frame(pca_med$scores)
pca_scores$ices_rect <- pop(pca_data)

# Load ICES area data
ices_tab <- read.csv("data/ices_area_rectangle.csv", sep = ",")
colnames(ices_tab) <- c("ices_rect", "ices_area")


# Merge PCA scores with ICES area information
pca_med_combined <- left_join(pca_scores, ices_tab, by = "ices_rect")

# Create Figure 5a (PCA plot including Mediterranean)
p1 <- ggplot(pca_med_combined, aes(x = PC1, y = PC2, colour = ices_area)) +
  geom_point(size = 2) +
  labs(x = paste0("PC1 (", PC1_variance, "% variance)"),
       y = paste0("PC2 (", PC2_variance, "% variance)")) +
  stat_ellipse(level = 0.95, size = 1)

# Save plot
pdf("output/pca/pca_med_cat.pdf")
print(p1)
dev.off()

# ---- PCA for Atlantic-only dataset (Figure 5b) ----

# Load Atlantic-only dataset
atl_data <- gl.load("data/bass_reduced_ices_spawning")

# Perform PCA
pca_atl <- glPca(atl_data, nf = 2)

# Compute variance explained
eig_total_atl <- sum(pca_atl$eig)
PC1_variance_atl <- formatC((pca_atl$eig[1] / eig_total_atl) * 100, format = "f", digits = 2)
PC2_variance_atl <- formatC((pca_atl$eig[2] / eig_total_atl) * 100, format = "f", digits = 2)

# Convert PCA scores to dataframe and add population labels
pca_atl_scores <- as.data.frame(pca_atl$scores)
pca_atl_scores$ices_rect <- pop(atl_data)

# Merge PCA scores with ICES area data
pca_atl_combined <- left_join(pca_atl_scores, ices_tab, by = "ices_rect")

# Create Figure 5b (PCA plot for Atlantic)
p2 <- ggplot(pca_atl_combined, aes(x = PC1, y = PC2, colour = ices_area)) +
  geom_point(size = 2) +
  labs(x = paste0("PC1 (", PC1_variance_atl, "% variance)"),
       y = paste0("PC2 (", PC2_variance_atl, "% variance)")) +
  stat_ellipse(level = 0.95, size = 1)

# Save plot
pdf("output/pca/pca_spawning_atl_ices_area.pdf")
print(p2)
dev.off()
