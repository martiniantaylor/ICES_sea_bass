
######################################################
##outflank analyses using dartR
#identifies outliers in atl & med and just atl pops using outflank 


#Martin Taylor (martin.taylor@uea.ac.uk)
######################################################

#load libraries
#devtools::install_github("whitlock/OutFLANK")
library(OutFLANK)
library(qvalue)
library(dartR)

#import file for atlantic pops (excludes portugal S)
outflank_atl <- readRDS("data/outflank/spawning_incl_glport") 
#med
outflank_med <- readRDS("data/outflank/spawning_incl_med") 

#check dartR compliance
outflank_med<-gl.compliance.check(outflank_med)
outflank_atl<-gl.compliance.check(outflank_atl)

#run outflank using dartR wrapper
out<-gl.outflank(
  outflank,
  plot = TRUE,
  LeftTrimFraction = 0.05,
  RightTrimFraction = 0.05,
  Hmin = 0.1,
  qthreshold = 0.05
)


out_med<-gl.outflank(
  outflank_med,
  plot = TRUE,
  LeftTrimFraction = 0.05,
  RightTrimFraction = 0.05,
  Hmin = 0.1,
  qthreshold = 0.05
)

#extract list of outliers snps from outfile
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

# Plot He versus Fst
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

#end of outFLANK
########################################################


ggplot(outliers_med,aes(x=-log10(qvalues),y=FST)) +
  geom_point(pch=21, size=2)+ 
#  scale_fill_manual(name="Selection",values=c("white","red","orange"))+ 
  labs(x="Log(q-value)")+ 
  labs(y="Fst")+ 
  labs(title="Outflank output spawning + Portugal N")+
  theme_classic()


#histogram of p values
hist(out_med$outflank$results$pvaluesRightTail)

med_out <- out_med$outflank$results
plot(med_out$He, med_out$FST, pch=19, col=rgb(0,0,0,0.1))
points(P1$He[med_out], P1$FST[med_out], col="blue")


#read_table output file

#outflank_res<-read.table("output/outflank/spawning_port.txt")

#colnames(outflank_res)<-c("LocusName",        "He",          "FST",            "T1",           "T2",  "FSTNoCorr", "T1NoCorr",     "T2NoCorr", "meanAlleleFreq", "indexOrder", "GoodH",   "qvalues","pvalues", "pvaluesRightTail", "OutlierFlag")

#plot for atlantic pops only
levels(out$results$OutlierFlag[1,])

outflak_df<-dplyr::select(outflank_res[1,2]) 

outflank_df<-as.data.frame(-log10(outflank_res$outflank.results.qvalues))
outflank_df$fst<-outflank_res$outflank.results.FST
colnames(outflank_df)<-c("qval","fst")

range(outflank_res$outflank.results.qvalues, na.rm=TRUE)


ggplot(outflank_df,aes(x=qval,y=fst)) +
  geom_point(pch=21, size=2)+ 
  scale_fill_manual(name="Selection",values=c("white","red","orange"))+ 
  labs(x="Log(q-value)")+ 
  labs(y="Fst")+ 
  labs(title="Outflank output spawning + Portugal N")+
  theme_classic()


#plot for atlantic + med  pops 
levels(out_med$outflank$results$OutlierFlag[1,])

outflank_df<-dplyr::select(outflank_res[1,2]) 

outflank_df<-as.data.frame(-log10(outflank_res$outflank.results.qvalues))
outflank_df$fst<-outflank_res$outflank.results.FST
colnames(outflank_df)<-c("qval","fst")

range(outflank_res$outflank.results.qvalues, na.rm=TRUE)


ggplot(outflank_df,aes(x=qval,y=fst)) +
  geom_point(pch=21, size=2)+ 
  scale_fill_manual(name="Selection",values=c("white","red","orange"))+ 
  labs(x="Log(q-value)")+ 
  labs(y="Fst")+ 
  labs(title="Outflank output spawning + Portugal N")+
  theme_classic()




######differences in snps identified by outflank and bayescan and pca adapt

#the bayescan object comes from the bayscan outlier script.
intesect1<-intersect(outliers_med2,outliers_med_bs_vec)

#18 loci in common
#[1] "AX-172301657_C" "AX-172281702_T" "AX-172283765_T" "AX-172318804_A"
#[5] "AX-172316821_G" "AX-172271956_G" "AX-172324113_T" "AX-172300868_G"
#[9] "AX-172276453_G" "AX-172296367_C" "AX-172308090_T" "AX-172320121_C"
#[13] "AX-172280467_T" "AX-172294467_A" "AX-172322750_A" "AX-172287166_T"
#[17] "AX-172296775_G" "AX-172274329_A"

intersect_outflank_pca_adapt<-intersect(outliers_med2,pcaasapt_outliers_vec)

length(intersect_outflank_pca_adapt)

intersect<-intersect(outliers_med_bs_vec,intersect_outflank_pca_adapt)

#import outliers
outflank_outliers<-readRDS("data/LEA_input/outlier_snps_outflank_red.gl_df.geno")

#all outliers from different methods
#bayescan snps are from V2 script
all_outliers<-c(outliers_med2,outliers_med_bs_vec,pcaasapt_outliers_vec)
length(all_outliers)

#get only unique ones
all_outliers_unique<-unique(all_outliers)

length(all_outliers_unique)

# Creating SNPs under selection object for all outlying snps------------------------------------
#this uses original genlight object used at start of pcadapt process
outlier_snps_all <- outflank_med[,all_outliers_unique]


#save to file as genlight object
gl.save(outlier_snps_all, "output/pca_adapt/sp_outlying_loci_all.gl", verbose = NULL)

