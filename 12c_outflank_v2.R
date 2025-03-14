######################################################
##outflank analyses using dartR
#identifies outliers in Atl & Med and just Atl pops using outflank 


#Martin Taylor (martin.taylor@uea.ac.uk)
######################################################

#load libraries
#devtools::install_github("whitlock/OutFLANK")
library(OutFLANK)
library(qvalue)
library(dartR)

#import file for atlantic pops (excludes portugal S)
outflank_atl <- gl.load("data/outflank/spawning_incl_glport.gl") 
#med
outflank_med <- gl.load("data/outflank/spawning_incl_med.gl") 

#check dartR compliance
outflank_med<-gl.compliance.check(outflank_med)
outflank_atl<-gl.compliance.check(outflank_atl)

#run outflank using dartR wrapper
#this is just Atlantic
out<-gl.outflank(
  outflank,
  plot = TRUE,
  LeftTrimFraction = 0.05,
  RightTrimFraction = 0.05,
  Hmin = 0.1,
  qthreshold = 0.05
)

#including MED
out_med<-gl.outflank(
  outflank_med,
  plot = TRUE,
  LeftTrimFraction = 0.05,
  RightTrimFraction = 0.05,
  Hmin = 0.1,
  qthreshold = 0.05
)

#extract list of outliers snps from outfile for MED
df_med<-out_med$outflank$results

# Print number of outliers (TRUE)
df_med$OutlierFlag %>% summary

#Mode   FALSE    TRUE 
#logical   39469    1553 

# Convert Fsts <0 to zero
df_med$FST[df_med$FST < 0] = 0 

# Italic labels
fstlab = expression(italic("F")[ST])
hetlab = expression(italic("H")[e])

# Plot He versus Fst for MED vs Atlantic
med_outliers<-ggplot(data = df_med)+
  geom_point(aes(x = He, y = FST, colour = OutlierFlag))+
  scale_colour_manual(values = c("black","red"), labels = c("Neutral SNP","Outlier SNP"))+
  ggtitle("OutFLANK outlier test")+
  xlab(hetlab)+
  ylab(fstlab)+
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 15, face = "bold")
  )

pdf("output/outflank/he_fst_med.pdf")
med_outliers
dev.off()


#filter for TRUE
outliers_med<-df_med[df_med$OutlierFlag == "TRUE", ]
#list of outlier loci
outliers_med2 <- out_med$outflank$results$LocusName[out_med$outflank$results$OutlierFlag == TRUE]
length(outliers_med2) #1553

#plot snps - full distribution
OutFLANKResultsPlotter(out_med$outflank, withOutliers = TRUE, NoCorr = TRUE, Hmin = 0.1, 
                       binwidth = 0.005, Zoom = FALSE, RightZoomFraction = 0.05, 
                       titletext = NULL)
#zoom on right tail
OutFLANKResultsPlotter(out_med$outflank, withOutliers = TRUE,
                       NoCorr = TRUE, Hmin = 0.1, binwidth = 0.001, Zoom =
                         TRUE, RightZoomFraction = 0.15, titletext = NULL)

#save genlight with just outliers
outlier_snps_outflank.gl <- gl.keep.loc(outflank_med, loc.list=outliers_med2)

#save outlier snps to file as genlight object
gl.save(outlier_snps_outflank.gl, "output/pca_adapt/sp_outlying_loci_outflank.gl", verbose = NULL)

#save locus names to file
outflank_outliers_vec<-outlier_snps_outflank.gl$loc.names

#####Atlantic outliers only
df_atl<-out$outflank$results

df_atl$OutlierFlag %>% summary

#Mode   FALSE    TRUE 
#logical   41021       1 

# Convert Fsts <0 to zero
df_atl$FST[df_atl$FST < 0] = 0 

# Italic labels
fstlab = expression(italic("F")[ST])
hetlab = expression(italic("H")[e])

# Plot He versus Fst
atl_outliers<-ggplot(data = df_atl)+
  geom_point(aes(x = He, y = FST, colour = OutlierFlag))+
  scale_colour_manual(values = c("black","red"), labels = c("Neutral SNP","Outlier SNP"))+
  ggtitle("OutFLANK outlier test - Atlantic")+
  xlab(hetlab)+
  ylab(fstlab)+
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 15, face = "bold")
  )

pdf("output/outflank/he_fst_atl.pdf")
atl_outliers
dev.off()

#filter for TRUE
outliers_atl<-df_atl[df_atl$OutlierFlag == "TRUE", ]
#list of outlier loci
outliers_atl2 <- out$outflank$results$LocusName[out$outflank$results$OutlierFlag == TRUE]
length(outliers_atl2) #1

#plot snps
OutFLANKResultsPlotter(out$outflank, withOutliers = TRUE, NoCorr = TRUE, Hmin = 0.1, 
                       binwidth = 0.005, Zoom = FALSE, RightZoomFraction = 0.05, 
                       titletext = NULL)

utils.outflank.plotter(out$outflank, withOutliers=T, NoCorr=TRUE, Hmin=0.1, binwidth=0.005, Zoom=FALSE, RightZoomFraction=0.05, titletext=NULL )

#
