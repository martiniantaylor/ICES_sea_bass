################################################################################
# Admixture / Distance re-analysis for Robinet Data
# 
# This script calculates the proportion of missing data for each individual,
# generates histograms of missing data, and plots the relationship between 
# Mediterranean ancestry and geographic distance from SINE for SNPs identified 
# in the Robinet (2000) dataset
#
# Author: Martin Taylor (martin.taylor@uea.ac.uk)
################################################################################

# Load necessary libraries
library(adegenet)  # For genetic data analysis
library(dplyr)     # For data manipulation
library(ggplot2)   # For data visualization

# Load SNP data
# Main dataset - 21 populations and 1012 loci
geno_827ATL_1012loci_21sreg <- read.genetix("data/robinet_data/FID_IID_DLAB1012loci_827ind_21reg.gtx")
#get names
"data/robinet_data/noms_827DLAB_ICESNAME_sreg_xy.csv" %>%
  read.csv() %>% as_tibble() -> noms_827DLAB

as.factor(noms_827DLAB$CLST) -> geno_827ATL_1012loci_21sreg$pop

# Calculate proportion of loci typed for each individual
prop_typed <- propTyped(geno_827ATL_1012loci_21sreg, by = c("ind"))
missing <- stack(prop_typed)
names(missing) <- c("typed", "specimen")

# Import SINE distances
distances <- read.csv("data/robinet_daat/dist_SINE_long_lat_31reg.csv")

# Generate histogram of missing data
robinet_histo <- ggplot(missing, aes(x = typed)) + 
  geom_histogram() +
  labs(title = "Histogram of Proportion of Loci Typed per Individual in Robinet Data", 
       x = "Proportion Typed", y = "Count") +
  theme_bw()

# Save histogram as PDF
pdf("histo_robinet_all.pdf", width = 6, height = 3)
print(robinet_histo)
dev.off()

# Load admixture data from Robinet paper (31 site data)
admix <- read.csv("data/robinet_data/res_admixture_FID_IID_1012loci_761ind_31CLST.csv")


# Merge missingness, admixture, and distance data
data_combined <- missing %>%
  left_join(admix, by = "specimen") %>%
  left_join(distances, by = "CLST") %>%
  mutate(CLST = factor(CLST, levels = c("SINE", "PENI", "PORT", "VIGO", "CORU", "ASTU1", "ASTU2", "CANT", 
                                        "AQUI1", "AQUI2", "CHAR4", "LOIR4", "LOIR5", "LOIR3", "PBRE2", 
                                        "PBRE1", "PBRE3", "CRNW3", "GONB1", "CRNW4", "CRNW2", "IRIS", 
                                        "BSEI1", "BSEI2", "ISWI2", "WALE2", "DOVE1", "SGEO", "DOVE3", 
                                        "ANGL", "NEUK3"))) %>%
  na.omit()  # Remove individuals with missing admixture scores

##make individual missingness plots

# Plot missing data vs MED ancestry
missing_plot_all <- ggplot(data_combined, aes(x = CLST, y = Q2, col = typed)) +
  geom_point(position = position_jitter(seed = 1, width = 0.2)) +
  labs(title = "Missing Data vs MED Ancestry in Robinet Data", 
       x = "Sample Site", y = "MED Admixture Proportion") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save plot as PDF
pdf("missingness_robinet_all.pdf", width = 6, height = 3)
print(missing_plot_all)
dev.off()

# Filter data to keep individuals genotyped > 95% markers
data_combined_95 <- data_combined %>% filter(typed >= 0.95)

# Plot filtered data vs MED ancestry
missing_plot_all_filt95 <- ggplot(data_combined_95, aes(x = CLST, y = Q2, col = typed)) +
  geom_point(position = position_jitter(seed = 1, width = 0.2)) +
  labs(title = "Missing Data vs MED Ancestry in Robinet Data (Filtered 95%)", 
       x = "Sample Site", y = "MED Admixture Proportion") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save plot as PDF
pdf("missingness_robinet_all_filt95.pdf", width = 6, height = 3)
print(missing_plot_all_filt95)
dev.off()

##-------MED admixture vs geographic distance plots---------

# Calculate mean and standard deviation for Q2 (filtered 95%)
q2_mean_all_95filt <- data_combined_95 %>%
  group_by(CLST) %>%
  summarise_at(vars(Q2), list(mean = mean, sd = sd)) %>%
  as.data.frame()

# Join distances with filtered MED ancestry
med_dist_95 <- left_join(distances, q2_mean_all_95filt)

# Add geographic regions for coloring
med_dist_95$geo <- c("North", "Biscay", "Biscay", "Biscay", "Biscay", "North", "North", "Biscay", 
                     "Biscay", "Portugal", "North", "North", "North", "North", "North", "North", 
                     "North", "North", "Biscay", "Biscay", "Biscay", "North", "Biscay", "Biscay", 
                     "Biscay", "Portugal", "Portugal", "North", "Portugal", "Portugal", "North")

# Plot distance to SINE vs MED ancestry (filtered 95%)
#this is the file with indivs with >5% missing data removed
robinet_95_redo <- ggplot(med_dist_95, aes(x = dist_SINE, y = mean, color = geo, label = CLST)) +
  geom_point() + 
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = .3) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "MED Ancestry in Robinet Data (Filtered 95%)", 
       x = "Distance to SINE (km)", y = "MED Admixture Proportion") +
  geom_text(hjust = 1, vjust = 2, angle = 90, size = 2) +
  theme_bw() +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# Save plot as PDF
pdf("robinet_all_MED_filt95.pdf", width = 6, height = 3)
print(robinet_95_redo)
dev.off()

# Summarize means and standard deviations for Q2 (unfiltered)
q2_mean_std_unfilt <- data_combined %>%
  group_by(CLST) %>%
  summarise_at(vars(Q2), list(mean = mean, sd = sd)) %>%
  as.data.frame()

# Join distances with unfiltered MED ancestry
med_dist_unfilt <- left_join(distances, q2_mean_std_unfilt)

# Add geographic regions for coloring
med_dist_unfilt$geo <- c("North", "Biscay", "Biscay", "Biscay", "Biscay", "North", "North", "Biscay", 
                         "Biscay", "Portugal", "North", "North", "North", "North", "North", "North", 
                         "North", "North", "Biscay", "Biscay", "Biscay", "North", "Biscay", "Biscay", 
                         "Biscay", "Portugal", "Portugal", "North", "Portugal", "Portugal", "North")

# Plot distance to SINE vs MED ancestry (unfiltered)
#this is the original data from robinet
robinet_original <- ggplot(med_dist_unfilt, aes(x = dist_SINE, y = mean, color = geo, label = CLST)) +
  geom_point() + 
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = .3) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "MED Ancestry in Robinet Data (Unfiltered)", 
       x = "Distance to SINE (km)", y = "MED Admixture Proportion") +
  geom_text(hjust = 1, vjust = 2, angle = 90, size = 2) +
  theme_bw() +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# Save plot as PDF
pdf("robinet_all_MED_unfiltered.pdf", width = 6, height = 3)
print(robinet_original)
dev.off()