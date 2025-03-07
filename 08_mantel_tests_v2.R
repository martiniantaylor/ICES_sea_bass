#checked and ready for pub
#######Mantel tests of feeding and spawning samples#######

# Load libraries
library(ade4)      # For mantel.rtest
library(ggplot2)   # For plotting
library(reshape2)  # For reshaping data

# Read Fst matrices - generated and saved in fst pairwise script,
Fst_matrix_spawn <- as.dist(read.delim("data/reordered_matrix_spawn.txt", header = TRUE, sep=" "))
Fst_matrix_feed <- as.dist(read.delim("data/reordered_matrix_feed_fst.txt", header = TRUE, sep=" "))

# Read geographic distance matrices - generated and saved in pairwise geographic distance script
dist_feed <- as.dist(read.delim("output/sea_distances/Contour_PairwiseMarineDistances_feeding.txt", sep = ",", header = TRUE))
dist_spawn <- as.dist(read.delim("output/sea_distances/Contour_Pairwise_Marine_Distances_spawning.txt", sep = ",", header = TRUE))

# Replace NA values with 0
dist_feed[is.na(dist_feed)] <- 0
dist_spawn[is.na(dist_spawn)] <- 0

# Log-transform geographic distances
dist_feed_log <- log(dist_feed)
dist_spawn_log <- log(dist_spawn)

# Run Mantel tests
mantel_feed <- mantel.rtest(Fst_matrix_feed, dist_feed_log, nrepet = 10000)
mantel_spawn <- mantel.rtest(Fst_matrix_spawn, dist_spawn_log, nrepet = 10000)

# Print Mantel test results
print(mantel_feed)
print(mantel_spawn)


# Function to prepare data for ggplot
prepare_plot_data <- function(dist_matrix, fst_matrix) {
  dist_df <- melt(dist_matrix)
  fst_df <- melt(fst_matrix)
  fst_fst2 <- fst_df$value / (1 - fst_df$value)  # Compute Fst/(1-Fst)
  data.frame(dist = dist_df$value, fst = fst_df$value, fst2 = fst_fst2)
}

# Prepare data for plotting
feeding_data <- prepare_plot_data(dist_feed_log, Fst_matrix_feed)
spawning_data <- prepare_plot_data(dist_spawn_log, Fst_matrix_spawn)


#ggplot for feeding data
feeding_plot<-ggplot(feeding_data, aes(x = dist, y = fst2)) + 
  geom_point() +
  stat_smooth(method = "lm", col = "black", se = FALSE, linetype = "dashed")+
  ggtitle(expression(atop("F"[ST]~"vs log sea distance feeding pops")))+
  labs( x = "log sea distance (km)", y = "Fst/1-Fst")+
  theme_bw()

#save ggplot in pdf
pdf(file="output/fst_analysis/genetic_distance_geographic_dist_atl_feeding.pdf",width = 8, # The width of the plot in inches
    height = 4) # The height of the plot in inches)))
feeding_plot
dev.off()

#ggplot for spawning data

spawning_plot<-ggplot(spawning_data, aes(x = dist, y = fst2)) + 
  geom_point() +
  stat_smooth(method = "lm", col = "black", se = FALSE, linetype = "dashed")+
  ggtitle(expression(atop("F"[ST]~"vs log sea distance spawning pops")))+
  labs( x = "log Sea distance (km)", y = "Fst/1-Fst")+
  theme_bw()

pdf(file="output/fst_analysis/genetic_distance_geographic_dist_atl_spawning.pdf",width = 8, # The width of the plot in inches
    height = 4) # The height of the plot in inches)))
spawning_plot
dev.off()
