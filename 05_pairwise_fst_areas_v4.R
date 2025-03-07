##script for pairwise fst calculation and associated heatmaps
##see earlier V3 of script for pairwise fst incl. MED /Portugal

# Load necessary libraries
library(dartR)
library(reshape)
library(ggplot2)
library(ggtext)
library(lmPerm)
library(dplyr)


#import data file
#spawning gl - 41022 loci, 286 indivs
bass_reduced_ices_spawning<-gl.load("data/bass_reduced_ices_spawning")

#should be atlantic pops only - no MED or port
bass_reduced_ices_spawning$pop

#use dartR to make fst matrix excl med 
dartr_fst_spawning<-gl.fst.pop(bass_reduced_ices_spawning, nboots=1000)

#save  fst matrices to file
saveRDS(dartr_fst_spawning, file="output/fst_analysis/spawning_fst.RData")
#read file
dartr_fst_spawning<-readRDS("output/fst_analysis/spawning_fst.RData")

#extract pairwise fsts
dartr_fst_pairwise_spawn<-dartr_fst_spawning$Fsts
#extract pairwise p-values 
dartr_pval_pairwise_spawn<-dartr_fst_spawning$Pvalues

#get list of spawning populations > 9 individuals
summary_spawn<-summary(bass_snps_spawning$pop)
pops_spawn<-summary_spawn[summary_spawn>9]
new_order_spawn<-names(pops_spawn)

#reorder FST
#give matrix size
n=15
# Initialize a matrix to store the reordered lower triangular part with NAs
reordered_matrix_spawn <- matrix(NA, nrow = n, ncol = n, dimnames = list(new_order_spawn, new_order_spawn))

# Reorder the lower triangular part while preserving NAs
for (i in 1:n) {
  for (j in 1:i) {
    if (!is.na(dartr_fst_pairwise_spawn[new_order_spawn[i], new_order_spawn[j]])) {
      reordered_matrix_spawn[i, j] <- dartr_fst_pairwise_spawn[new_order_spawn[i], new_order_spawn[j]]
    } else {
      reordered_matrix_spawn[i, j] <- dartr_fst_pairwise_spawn[new_order_spawn[j], new_order_spawn[i]]
    }
  }
}

#save matrix to file for later on
write.table(reordered_matrix_spawn, file="data/reordered_matrix_spawn.txt", row.names=TRUE, col.names=TRUE)

#reorder p-values
new_order_spawn_pvals <- pops_keep_spawn
#give matrix size
n=15
# Initialize a matrix to store the reordered lower triangular part with NAs
reordered_matrix_spawn_pvals <- matrix(NA, nrow = n, ncol = n, dimnames = list(new_order_spawn_pvals, new_order_spawn_pvals))


# Reorder the lower triangular part while preserving NAs
for (i in 1:n) {
  for (j in 1:i) {
    if (!is.na(dartr_pval_pairwise_spawn[new_order_spawn_pvals[i], new_order_spawn_pvals[j]])) {
      reordered_matrix_spawn_pvals[i, j] <- dartr_pval_pairwise_spawn[new_order_spawn_pvals[i], new_order_spawn_pvals[j]]
    } else {
      reordered_matrix_spawn_pvals[i, j] <- dartr_pval_pairwise_spawn[new_order_spawn_pvals[j], new_order_spawn_pvals[i]]
    }
  }
}

#save matric to file for later on
write.table(reordered_matrix_spawn_pvals, file="data/reordered_matrix_spawn_pvals.txt", row.names=TRUE, col.names=TRUE)


#convert fsts to data frame- can use melt or cor_gather
fst.mat_test_sp<-reshape::melt(reordered_matrix_spawn)
str(fst.mat_test_sp)
fst.mat_test_sp$value<-as.numeric(fst.mat_test_sp$value)

#convert pvalues to data frame- can use melt or cor_gather
fst.mat_test_sp_pvals<-reshape::melt(reordered_matrix_spawn_pvals)
str(fst.mat_test_sp_pvals)
fst.mat_test_sp_pvals$value<-as.numeric(fst.mat_test_sp_pvals$value)
tests<-(nrow(fst.mat_test_sp_pvals)-sum(is.na(fst.mat_test_sp_pvals$value)))

fst.mat_test_sp_pvals$adj<-p.adjust(fst.mat_test_sp_pvals$value,method="bonferroni", n=tests)

fst.mat_test_sp_pvals$stars <- cut(fst.mat_test_sp_pvals$adj, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))  # Create column of significance labels

# Convert minus values to zero
fst.mat_test_sp$value[fst.mat_test_sp$value < 0] = 0

# Fst italic label
fst.label = expression(italic("F")[ST])

#plot heatmap spawning season. ices rectangles.
library(ggtext)
fst_spawning<-ggplot(data = fst.mat_test_sp, aes(x = X2, y = X1, fill = value))+
  geom_tile(color = "black")+ scale_fill_gradient(low = "white", high = "red", name=fst.label, na.value = "white") + 
  geom_text(label=fst.mat_test_sp_pvals$stars, color="black", size=3, vjust = -0.5) + 
  geom_text(label=round(fst.mat_test_sp$value, 4), color="black", size=2.5) + 
  ggtitle("Pairwise F<sub>ST</sub> spawning pops, WC (1984)")+
  labs( x = "ICES rectangle", y = "ICES rectangle") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 9, hjust = 1),
        axis.text.y = element_text(size = 9),
        plot.title = element_markdown(size = 12, lineheight = 1.2)) + 
  coord_fixed()

####for feeding populations

#import feeding data set gl - 41022 loci, 558 indivs
bass_reduced_ices_feeding<-gl.load("data/bass_reduced_ices_feeding")


bass_reduced_ices_feeding$pop

#use dartR to make fst matrix excl med 
dartr_fst_feeding<-gl.fst.pop(bass_reduced_ices_feeding, nboots=1000)

#save  fst matrices to file
saveRDS(dartr_fst_feeding, file="output/fst_analysis/feeding_fst.RData")
dartr_fst_feeding<-readRDS("output/fst_analysis/feeding_fst.RData")

dartr_fst_pairwise_feed<-dartr_fst_feeding$Fsts
dartr_pval_pairwise_feed<-dartr_fst_feeding$Pvalues

#this file comes from the file making script - for feeding pops
#first get list of feeding populations > 9 individuals
summary_feed<-summary(bass_snps_feeding$pop)
pops_feed<-summary_feed[summary_feed>9]
pops_keep_feed<-names(pops_feed)

#reorder the matrices into sensible order

#give matrix size
n=20
# Initialize a matrix to store the reordered lower triangular part with NAs
reordered_matrix_feed <- matrix(NA, nrow = n, ncol = n, dimnames = list(new_order_feed, new_order_feed))

for (i in 1:n) {
  for (j in 1:i) {
    if (!is.na(dartr_fst_pairwise_feed[new_order_feed[i], new_order_feed[j]])) {
      reordered_matrix_feed[i, j] <- dartr_fst_pairwise_feed[new_order_feed[i], new_order_feed[j]]
    } else {
      reordered_matrix_feed[i, j] <- dartr_fst_pairwise_feed[new_order_feed[j], new_order_feed[i]]
    }
  }
}

#save matrix to file for later on
write.table(reordered_matrix_feed, file="data/reordered_matrix_feed_fst.txt", row.names=TRUE, col.names=TRUE)


new_order_feed_pvals <- pops_keep_feed
#give matrix size
n=20
# Initialize a matrix to store the reordered lower triangular part with NAs
reordered_matrix_feed_pvals <- matrix(NA, nrow = n, ncol = n, dimnames = list(new_order_feed_pvals, new_order_feed_pvals))


# Reorder the lower triangular part while preserving NAs
for (i in 1:n) {
  for (j in 1:i) {
    if (!is.na(dartr_pval_pairwise_feed[new_order_feed_pvals[i], new_order_feed_pvals[j]])) {
      reordered_matrix_feed_pvals[i, j] <- dartr_pval_pairwise_feed[new_order_feed_pvals[i], new_order_feed_pvals[j]]
    } else {
      reordered_matrix_feed_pvals[i, j] <- dartr_pval_pairwise_feed[new_order_feed_pvals[j], new_order_feed_pvals[i]]
    }
  }
}

#save matrix to file for later on
write.table(reordered_matrix_feed_pvals, file="data/reordered_matrix_feed_pval.txt", row.names=TRUE, col.names=TRUE)


#convert to data frame- can use melt or cor_gather
fst.mat_test_f<-reshape::melt(reordered_matrix_feed)
str(fst.mat_test_f)
fst.mat_test_f$value<-as.numeric(fst.mat_test_f$value)

#convert pvalues to data frame- can use melt or cor_gather
fst.mat_test_feed_pvals<-reshape::melt(reordered_matrix_feed_pvals)
str(fst.mat_test_feed_pvals)
fst.mat_test_feed_pvals$value<-as.numeric(fst.mat_test_feed_pvals$value)
tests<-(nrow(fst.mat_test_feed_pvals)-sum(is.na(fst.mat_test_feed_pvals$value)))

fst.mat_test_feed_pvals$adj<-p.adjust(fst.mat_test_feed_pvals$value,method="bonferroni", n=tests)

fst.mat_test_feed_pvals$stars <- cut(fst.mat_test_feed_pvals$adj, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))  # Create column of significance labels




# Convert minus values to zero
fst.mat_test_f$value[fst.mat_test_f$value < 0] = 0

# Fst italic label
fst.label = expression(italic("F")[ST])

#plot heatmap spawning season. ices rectangles.

fst_feeding<-ggplot(data = fst.mat_test_f, aes(x = X2, y = X1, fill = value))+
  geom_tile(color = "black")+ scale_fill_gradient(low = "white", high = "red", name=fst.label, na.value = "white") + 
  geom_text(label=fst.mat_test_feed_pvals$stars, color="black", size=3,  vjust = -0.3) + 
  geom_text(label=round(fst.mat_test_f$value, 4), color="black", size=2) + 
  ggtitle(expression(atop("Pairwise FST feeding pops, WC (1984) Dartr")))+
  labs( x = "ICES rectangle", y = "ICES rectangle") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 9, hjust = 1),axis.text.y = element_text(size = 9)) + 
  coord_fixed()

pdf("output/fst_analysis/fst_heatmaps.pdf")
par(mfrow = c(1, 2))
fst_feeding
fst_spawning
dev.off()

#test fst between feeding and spawning--------------

#rbind two fst tibbles

#insert grouping variable
spawning<-fst.mat_test_sp
spawning$group<-"spawning"

feeding<-fst.mat_test_f
feeding$group<-"feeding"


test_fst<-rbind(spawning,feeding)

#removes rows with nas and missing data
fst_nonas<-test_fst[complete.cases(test_fst), ]

fst_nonas$group<-as.factor(fst_nonas$group)

#box plot of outputs
perm_test<-ggplot(fst_nonas, aes(x=group, y=value))+
  geom_boxplot()+
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.3)+
  ggtitle(expression(atop("Pairwise FST feeding samples vs Spawning samples")))+
  labs( x = "Season", y = "Pairwise FST") +
  theme_bw()


pdf(file="output/fst_analysis/fst_boxplot.pdf")
perm_test
dev.off()


#permutation test using lmperm
library(lmPerm)
summary(lmp(value~group,data=fst_nonas))

#[1] "Settings:  unique SS "
#
#Call:
#  lmp(formula = value ~ group, data = fst_nonas)
#
#Residuals:
#  Min         1Q     Median         3Q        Max 
#-4.218e-04 -2.671e-04 -7.772e-05  1.684e-04  1.486e-03 
#
#Coefficients:
#  Estimate Iter Pr(Prob)  
#group1 -4.982e-05 5000   0.0106 *
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 0.0003359 on 293 degrees of freedom
#Multiple R-Squared: 0.0199,	Adjusted R-squared: 0.01655 
#F-statistic: 5.949 on 1 and 293 DF,  p-value: 0.01532 



#Pairwise Fst for ices areas---------------------------------------------------------

#get ices area for each ices rectangle (prev saved in script 03)
ices_tab_uniq<-read.csv("data/ices_area_rectangle.csv", sep=",")

#get existing pops spawning
spawn_ices_rects<-data.frame(bass_reduced_ices_spawning$pop, stringsAsFactors = FALSE)
colnames(spawn_ices_rects)<-"ices_rect"
#left join spawning pops with ful list
spawning_areas<-left_join(spawn_ices_rects,ices_tab_uniq )
#save as file with no colnames or row names
write.table(spawning_areas, file="data/temp_files/spawning_areas_rename.txt", row.names = F, col.names = F, quote = FALSE, sep =",")

#rename populations using dartr function
bass_reduced_ices_area_spawning<-gl.recode.pop(bass_reduced_ices_spawning, pop.recode="data/temp_files/spawning_areas_rename.txt", recalc = TRUE, mono.rm = FALSE, verbose = 5)
#sanity check 5 - make sure that this ices area and not rectangles still
bass_reduced_ices_area_spawning$pop

#use dartR to make fst matrix excl med
dartr_fst_spawning_ices_areas<-gl.fst.pop(bass_reduced_ices_area_spawning, nboots=1000)

#extract pairwise fsts
dartr_fst_pairwise_spawn_areas<-dartr_fst_spawning_ices_areas$Fsts
#p values
dartr_pvalues_pairwise_spawn_areas<-dartr_fst_spawning_ices_areas$Pvalues


#get icea area for each pop
areas_keep_spawn<-as.vector(unique(bass_reduced_ices_area_spawning$pop))

#give matrix size
n=6
# Initialize a matrix to store the reordered lower triangular part with NAs
reordered_matrix_spawn_area <- matrix(NA, nrow = n, ncol = n, dimnames = list(areas_keep_spawn, areas_keep_spawn))


# Reorder the lower triangular part while preserving NAs
for (i in 1:n) {
  for (j in 1:i) {
    if (!is.na(dartr_fst_pairwise_spawn_areas[areas_keep_spawn[i], areas_keep_spawn[j]])) {
      reordered_matrix_spawn_area[i, j] <-dartr_fst_pairwise_spawn_areas[areas_keep_spawn[i], areas_keep_spawn[j]]
    } else {
      reordered_matrix_spawn_area[i, j] <- dartr_fst_pairwise_spawn_areas[areas_keep_spawn[j], areas_keep_spawn[i]]
    }
  }
}

#give matrix size
n=6
# Initialize a matrix to store the reordered lower triangular part with NAs
reordered_matrix_spawn_area_Pval <- matrix(NA, nrow = n, ncol = n, dimnames = list(areas_keep_spawn, areas_keep_spawn))


# Reorder the lower triangular part while preserving NAs
for (i in 1:n) {
  for (j in 1:i) {
    if (!is.na(dartr_pvalues_pairwise_spawn_areas[areas_keep_spawn[i], areas_keep_spawn[j]])) {
      reordered_matrix_spawn_area_Pval[i, j] <-dartr_pvalues_pairwise_spawn_areas[areas_keep_spawn[i], areas_keep_spawn[j]]
    } else {
      reordered_matrix_spawn_area_Pval[i, j] <- dartr_pvalues_pairwise_spawn_areas[areas_keep_spawn[j], areas_keep_spawn[i]]
    }
  }
}





#convert to data frame- using melt
fst.mat_test_sp_area<-reshape::melt(reordered_matrix_spawn_area)
str(fst.mat_test_sp_area)
fst.mat_test_sp_area$value<-as.numeric(fst.mat_test_sp_area$value)

#convert pvalues to data frame- can use melt or cor_gather
fst.mat_test_sp_pvals<-reshape::melt(reordered_matrix_spawn_area_Pval)
str(fst.mat_test_sp_pvals)
fst.mat_test_sp_pvals$value<-as.numeric(fst.mat_test_sp_pvals$value)
tests<-(nrow(fst.mat_test_sp_pvals)-sum(is.na(fst.mat_test_sp_pvals$value)))

fst.mat_test_sp_pvals$adj<-p.adjust(fst.mat_test_sp_pvals$value,method="bonferroni", n=tests)

fst.mat_test_sp_pvals$stars <- cut(fst.mat_test_sp_pvals$adj, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))  # Create column of significance labels



# Convert minus values to zero
fst.mat_test_sp_area$value[fst.mat_test_sp_area$value < 0] = 0


# Fst italic label
fst.label = expression(italic("F")[ST])

#plot heatmap spawning season. ices areas. added theme parts to remove grey border to heatmap
options(scipen = 999)
fst_spawning_area_ices<-ggplot(data = fst.mat_test_sp_area, aes(x = X2, y = X1, fill = value))+
  geom_tile(color = "black")+ scale_fill_gradient(low = "white", high = "red", name=fst.label, na.value = "white") + 
  geom_text(label=fst.mat_test_sp_pvals$stars, color="black", size=3, vjust = -0.7) + 
  geom_text(label=round(fst.mat_test_sp_area$value, 4), color="black", size=3) + 
  ggtitle(expression(atop("Pairwise FST spawning pops, WC (1984) Dartr")))+
  labs( x = "ICES division", y = "ICES division") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 9, hjust = 1),axis.text.y = element_text(size = 9), panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    coord_fixed()

pdf(file="output/fst_analysis/fst_spawning_area_ices_sep_2024.pdf")
fst_spawning_area_ices
dev.off()

#fst heatmap ices areas/divisions for feeding populations

#get existing pops feeding
feed_ices_rects<-data.frame(bass_reduced_ices_feeding$pop, stringsAsFactors = FALSE)
colnames(feed_ices_rects)<-"ices_rect"
#left join spawning pops with ful list
feeding_areas<-left_join(feed_ices_rects,ices_tab_uniq )
#save as file with no colnames or row names
write.table(feeding_areas, file="data/temp_files/feeding_areas_rename.txt", row.names = F, col.names = F, quote = FALSE, sep =",")

#rename populations using dartr function
bass_reduced_ices_area_feeding<-gl.recode.pop(bass_reduced_ices_feeding, pop.recode="data/temp_files/feeding_areas_rename.txt", recalc = TRUE, mono.rm = FALSE, verbose = 5)
#sanity check 5 - make sure that this iceas area and not rectangles still
bass_reduced_ices_area_feeding$pop

#use dartR to make fst matrix excl med
dartr_fst_feeding_ices_areas<-gl.fst.pop(bass_reduced_ices_area_feeding, nboots=1000)

#extract pairwise fsts
dartr_fst_pairwise_feed_areas<-dartr_fst_feeding_ices_areas$Fsts
#p values
dartr_pvalues_pairwise_feed_areas<-dartr_fst_feeding_ices_areas$Pvalues


#this file comes from the file making script - for spawning pops
areas_keep_feed<-as.vector(unique(bass_reduced_ices_area_feeding$pop))

#give matrix size
n=9
# Initialize a matrix to store the reordered lower triangular part with NAs
reordered_matrix_feed_area <- matrix(NA, nrow = n, ncol = n, dimnames = list(areas_keep_feed, areas_keep_feed))


# Reorder the lower triangular part while preserving NAs
for (i in 1:n) {
  for (j in 1:i) {
    if (!is.na(dartr_fst_pairwise_feed_areas[areas_keep_feed[i], areas_keep_feed[j]])) {
      reordered_matrix_feed_area[i, j] <-dartr_fst_pairwise_feed_areas[areas_keep_feed[i], areas_keep_feed[j]]
    } else {
      reordered_matrix_feed_area[i, j] <- dartr_fst_pairwise_feed_areas[areas_keep_feed[j], areas_keep_feed[i]]
    }
  }
}

#give matrix size
n=9
# Initialize a matrix to store the reordered lower triangular part with NAs
reordered_matrix_feed_area_Pval <- matrix(NA, nrow = n, ncol = n, dimnames = list(areas_keep_feed, areas_keep_feed))


# Reorder the lower triangular part while preserving NAs
for (i in 1:n) {
  for (j in 1:i) {
    if (!is.na(dartr_pvalues_pairwise_feed_areas[areas_keep_feed[i], areas_keep_feed[j]])) {
      reordered_matrix_feed_area_Pval[i, j] <-dartr_pvalues_pairwise_feed_areas[areas_keep_feed[i], areas_keep_feed[j]]
    } else {
      reordered_matrix_feed_area_Pval[i, j] <- dartr_pvalues_pairwise_feed_areas[areas_keep_feed[j], areas_keep_feed[i]]
    }
  }
}





#convert to data frame- can use melt or cor_gather
fst.mat_test_feed_area<-reshape::melt(reordered_matrix_feed_area)
str(fst.mat_test_feed_area)
fst.mat_test_feed_area$value<-as.numeric(fst.mat_test_feed_area$value)

#convert pvalues to data frame- can use melt or cor_gather
fst.mat_test_feed_pvals<-reshape::melt(reordered_matrix_feed_area_Pval)
str(fst.mat_test_feed_pvals)
fst.mat_test_feed_pvals$value<-as.numeric(fst.mat_test_feed_pvals$value)
tests<-(nrow(fst.mat_test_feed_pvals)-sum(is.na(fst.mat_test_feed_pvals$value)))

fst.mat_test_feed_pvals$adj<-p.adjust(fst.mat_test_feed_pvals$value,method="bonferroni", n=tests)

fst.mat_test_feed_pvals$stars <- cut(fst.mat_test_feed_pvals$adj, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))  # Create column of significance labels



# Convert minus values to zero
fst.mat_test_feed_area$value[fst.mat_test_feed_area$value < 0] = 0


# Fst italic label
fst.label = expression(italic("F")[ST])

#plot heatmap spawning season. ices areas. added theme parts to remove grey border to heatmap
options(scipen = 999)
fst_feeding_area_ices<-ggplot(data = fst.mat_test_feed_area, aes(x = X2, y = X1, fill = value))+
  geom_tile(color = "black")+ scale_fill_gradient(low = "white", high = "red", name=fst.label, na.value = "white") + 
  geom_text(label=fst.mat_test_feed_pvals$stars, color="black", size=3, vjust = -0.7) + 
  geom_text(label=round(fst.mat_test_feed_area$value, 4), color="black", size=3) + 
  ggtitle(expression(atop("Pairwise FST feeding pops, WC (1984) Dartr")))+
  labs( x = "ICES division", y = "ICES division") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 9, hjust = 1),axis.text.y = element_text(size = 9), panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  coord_fixed()

pdf(file="output/fst_analysis/fst_feeding_area_ices_sep_2024.pdf")
fst_feeding_area_ices
dev.off()