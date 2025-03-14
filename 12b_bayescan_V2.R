########################################################################
# BayeScan Outlier Detection 

# This script runs BayeScan for outlier detection using SNP data from two datasets:
# - Atlantic spawning populations (including Portugal N)
# - Mediterranean populations
#
#Martin Taylor (martin.taylor@uea.ac.uk)
######################################################

# Load required libraries
library(LEA)       # Genomic data handling
library(dartR)     # SNP data processing
library(coda)      # MCMC diagnostics
library(ggplot2)   # Data visualization

# -----------------------------------------------------------------------------
# 1. Load SNP Data
# -----------------------------------------------------------------------------
# Load genlight objects for Atlantic and Mediterranean populations
for_bayescan_atl <- gl.load("data/outflank/spawning_incl_glport.gl")
for_bayescan_med <- gl.load("data/outflank/spawning_incl_med.gl")

# Check populations
print(levels(for_bayescan_atl$pop))
print(levels(for_bayescan_med$pop))

# Check dartR compliance
for_bayescan_atl <- gl.compliance.check(for_bayescan_atl)
for_bayescan_med <- gl.compliance.check(for_bayescan_med)

# -----------------------------------------------------------------------------
# 2. Convert GenLight Object to BayeScan Format
# -----------------------------------------------------------------------------
gl2bayescan(for_bayescan_atl, outfile = "bayescan_sp_atl.txt", outpath = "output/bayescan")
gl2bayescan(for_bayescan_med, outfile = "bayescan_sp_med.txt", outpath = "output/bayescan")

# -----------------------------------------------------------------------------
# 3. Run BayeScan (External Command Execution)
# -----------------------------------------------------------------------------
# Execute the following commands in the terminal:
# BayeScan settings: 
# - 5000 iterations, 10 thinning, 20 pilot runs, 50,000 burn-in, prior odds 100

# bayescan bayescan_sp_atl.txt -o bayescanLongRunOD100 n 5000 -thin 10 -nbp 20 -pilot 5000 -burn 50000 -pr_odds 100
# bayescan bayescan_sp_med.txt -o bayescanLongRun_medOD100 n 5000 -thin 10 -nbp 20 -pilot 5000 -burn 50000 -pr_odds 100

# -----------------------------------------------------------------------------
# 4. Evaluate Convergence (MCMC Diagnostics)
# -----------------------------------------------------------------------------
# Load MCMC output for Atlantic and Mediterranean runs
chain_atl <- read.table("output/bayescan/bayescanLongRunOD100.sel", header = TRUE)
chain_med <- read.table("output/bayescan/bayescanLongRun_medOD100.sel", header = TRUE)

# Source BayeScan plotting script (make sure plot_R.r exists in your directory)
source("scripts/plot_R.r")

# BayeScan plot FST distribution
plot_bayescan("output/bayescan/bayescanLongRunOD100_fst.txt", FDR = 0.05)
plot_bayescan("output/bayescan/bayescanLongRun_medOD100_fst.txt", FDR = 0.05)

# make  MCMC diagnostics function
mcmc_diag <- function(chain_data) {
  chain_data <- chain_data[, -1] # Remove first column (ID)
  chain_mcmc <- mcmc(chain_data, thin = 10)
  
  plot(chain_mcmc)                  # Visual check of chain mixing
  print(summary(chain_mcmc))         # Summary statistics
  print(autocorr.diag(chain_mcmc))   # Autocorrelation check
  print(effectiveSize(chain_mcmc))   # Effective sample size
  print(geweke.diag(chain_mcmc))     # Geweke diagnostic
  print(heidel.diag(chain_mcmc))     # Heidelberger and Welch diagnostic
}

par(mar=c(1,1,1,1))
# Run diagnostics for each dataset
mcmc_diag(chain_atl)
mcmc_diag(chain_med)

# -----------------------------------------------------------------------------
# 5. Analyze & Visualize Outliers
# -----------------------------------------------------------------------------

# Load FST results for Atlantic and Mediterranean
SNP_atl <- read.table("output/bayescan/bayescanLongRunOD100_fst.txt", header = TRUE)
SNP_med <- read.table("output/bayescan/bayescanLongRun_medOD100_fst.txt", header = TRUE)

# Assign SNP names
SNP_atl$snp <- for_bayescan_atl$loc.names
SNP_med$snp <- for_bayescan_med$loc.names

# make function to ensure q-values are numeric & round values for clarity
process_SNP_data <- function(SNP_data) {
  SNP_data$qval <- as.numeric(SNP_data$qval)
  SNP_data[SNP_data$qval <= 0.0001, "qval"] <- 0.0001 # Set minimum threshold
  #round values
  SNP_data$log10_PO <- round(SNP_data$log10.PO., 4) 
  SNP_data$qval <- round(SNP_data$qval, 4)
  SNP_data$alpha <- round(SNP_data$alpha, 4)
  SNP_data$fst <- round(SNP_data$fst, 6)
  
  # Classify selection type based on q-value and alpha
  SNP_data$SELECTION <- ifelse(SNP_data$alpha >= 0 & SNP_data$qval <= 0.05, "diversifying",
                               ifelse(SNP_data$alpha >= 0 & SNP_data$qval > 0.05, "neutral", "balancing"))
  
  SNP_data$SELECTION <- factor(SNP_data$SELECTION)
  return(SNP_data)
}

SNP_atl <- process_SNP_data(SNP_atl)
SNP_med <- process_SNP_data(SNP_med)

# Summary of selection types
print(xtabs(data = SNP_atl, ~SELECTION))
print(xtabs(data = SNP_med, ~SELECTION))

# Plot BayeScan results for each dataset

#ATL + MED
ggplot(SNP_med,aes(x=log10_PO,y=fst)) +
  geom_point(aes(fill=SELECTION), pch=21, size=2)+ 
  scale_fill_manual(name="Selection",values=c("white","red","orange"))+ 
  labs(x="Log(q-value)")+ 
  labs(y="Fst")+ 
  labs(title="Bayescan output spawning + Med")+
  theme_classic()

#ATL
ggplot(SNP_atl,aes(x=log10_PO,y=fst)) +
  geom_point(aes(fill=SELECTION), pch=21, size=2)+ 
  scale_fill_manual(name="Selection",values=c("white","red","orange"))+ 
  labs(x="Log(q-value)")+ 
  labs(y="Fst")+ 
  labs(title="Bayescan output spawning + Portugal N")+
  theme_classic()

# -----------------------------------------------------------------------------
# 6. Save Outliers for Further Analysis
# -----------------------------------------------------------------------------

# Extract outlier SNPs (diversifying selection) from Mediterranean dataset
outliers_med_bs <- SNP_med[SNP_med$SELECTION == "diversifying", ]
outliers_med_bs_vec <- outliers_med_bs$snp

# Filter genlight object to retain only outlier loci
outlier_snps_bayescan_gl <- gl.keep.loc(for_bayescan_med, loc.list = outliers_med_bs_vec)

# Save outlier SNPs as a genlight object
gl.save(outlier_snps_bayescan_gl, "output/bayescan/outlying_loci_bayescan.gl")

