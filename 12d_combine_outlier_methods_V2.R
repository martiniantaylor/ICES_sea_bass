################################################################################
# Combining SNPs Identified by Outflank, Bayescan, and PCAadapt
# 
# This script loads SNPs identified as outliers from three different selection 
# detection methods (Outflank, Bayescan, and pcadapt), merges them into a 
# single list, and finds SNPs common across all methods.
#
# Author: Martin Taylor (martin.taylor@uea.ac.uk)
################################################################################

# Load required libraries
library(dartR)  # For handling genlight objects

# ------------------------------------------------------------------------------
# 1. Load Outlier SNPs from Different Methods
# ------------------------------------------------------------------------------

# Load outliers detected by Bayescan
bayescan_outliers <- gl.load("output/bayescan/outlying_loci_bayescan.gl")
bayescan_outliers_vec <- bayescan_outliers$loc.names  # Extract SNP names

# Load outliers detected by Outflank
outflank_outliers <- gl.load("output/outflank/sp_outlying_loci_outflank.gl")
outflank_outliers_vec <- outflank_outliers$loc.names  # Extract SNP names

# Load outliers detected by PCAadapt
pcaadapt_outliers <- gl.load("output/pca_adapt/sp_outlying_loci_pca_adapt.gl")
pcaadapt_outliers_vec <- pcaadapt_outliers$loc.names  # Extract SNP names

# ------------------------------------------------------------------------------
# 2. Load the Full Dataset (for extracting selected SNPs later)
# ------------------------------------------------------------------------------

for_outliers_atl_med.gl <- gl.load("data/outflank/spawning_incl_med.gl") 

# ------------------------------------------------------------------------------
# 3. Merge All Outlier SNPs from Different Methods
# ------------------------------------------------------------------------------

# Combine all unique outlier SNPs detected by any method
all_outliers <- unique(c(pcaadapt_outliers_vec, outflank_outliers_vec, bayescan_outliers_vec))
cat("Total unique outlier SNPs:", length(all_outliers), "\n")

# Extract the identified outlier SNPs from the full dataset
outlier_snps_all <- for_outliers_atl_med.gl[, all_outliers]

# Save the combined outlier SNPs as a genlight object for future use and plotting
gl.save(outlier_snps_all, "output/outflank/sp_outlying_loci_all.gl", verbose = NULL)

# ------------------------------------------------------------------------------
# 4. Identify SNPs Common to All Three Methods
# ------------------------------------------------------------------------------

# Find SNPs that were detected as outliers by all three methods
common_snps <- Reduce(intersect, list(bayescan_outliers_vec, outflank_outliers_vec, pcaadapt_outliers_vec))

# Print the number of common SNPs
cat("Number of SNPs identified as outliers by all methods:", length(common_snps), "\n")

# Display the common SNPs
print(common_snps)

