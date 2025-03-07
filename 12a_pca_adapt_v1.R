######################################################
## PCAadapt Outlier identification on spawning pops
#saves pca dapt identified outliers as new genlight file


#Martin Taylor (martin.taylor@uea.ac.uk)
######################################################

#load libraries
library(LEA)         # Population genetics analysis
library(pcadapt)     # PCA-based adaptation
library(qvalue)      # Multiple hypothesis testing correction
library(dplyr)       # Data manipulation


#including med - this is 41022 locus data set incl portN and MED
pca_adapt_med <- readRDS("data/outflank/spawning_incl_med") 

#check pops
pca_adapt_med$pop

#convert to dataframe
pca_adapt_med.df <- as.data.frame(pca_adapt_med) 

# remove nas and replace with -9
pca_adapt_med.df[is.na(pca_adapt_med.df)] <- 9

#write lfmm file (uses LEA function)
write.lfmm(pca_adapt_med.df, "bass_snps_pca_adapt_med.lfmm")

#convert lfmm to pcadapt format

bs_pcaadapt <- read.pcadapt("bass_snps_pca_adapt_med.lfmm",
                            type = "lfmm")

#creating PCA in pcaadpt, using 10 PCs - as that captures the majority of variability
bs_pcaadapt_pca <- pcadapt(input = bs_pcaadapt, K=10)

#scree plot - shows a very long a shallow decline after steep cline from 1 to 2 PCs
plot(bs_pcaadapt_pca,option="screeplot")

res <- pcadapt(input = bs_pcaadapt, K = 2)
plot(res)

plot(res, option = "qqplot")

#shows that there are some area of the genome that have low recombination....
par(mfrow = c(2, 2))
for (i in 1:4)
  plot(res$loadings[, i], pch = 19, cex = .3, ylab = paste0("Loadings PC", i))

#split colours into atlantic, med and portugal
poplist.names <- c(rep("ATL", 286),rep("MED", 6),rep("PORT", 25))

#pca plot showing basins
plot(res,option="scores", pop = poplist.names)

#identifying statistical outliers using q
qval <- qvalue(bs_pcaadapt_pca$pvalues)$qvalues
#save loci that have a false discovery rate of less than 1%
outliers <- which(qval<0.01)
#how many outliers?
length(outliers)

#run this if error for figure margins being too large
par(mar=c(1,1,1,1))
#make histogram
hist(res$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange")

# Creating SNPs under selection object ------------------------------------
#this filters original genlight object used at start of pcadapt process
outlier_snps_pca_adapt <- pca_adapt_med[,outliers]

#save outlier snps to file as genlight object
gl.save(outlier_snps_pca_adapt, "output/pca_adapt/sp_outlying_loci_pca_adapt.gl", verbose = NULL)

#save locus names to file
pcaasapt_outliers_vec<-outlier_snps_pca_adapt$loc.names
