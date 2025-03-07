######################################################
#Pairwise population fst script using dartR
#Figure 1 in sea bass manuscript
#This will make combined heatmap using both feeding 
#and spawning data 
#see script ?? for pairwise ICES area heatmaps and 
#statistical test of differences among feeding and spawning 

#Martin Taylor (martin.taylor@uea.ac.uk)
######################################################

#use genlight file saved in script 05c - pop specific fst
#import feeding and spawning data set gl - 41022 loci, 901 indivs
bass_all_incl_med<-gl.load("data/combined_incl_med.gl")

#compliance check for dartr
bass_all_incl_med<-gl.compliance.check(bass_all_incl_med, verbose = 5)

#check pops
levels(bass_all_incl_med$pop)

#remove med and Portugal N and S pops
bass_all_excl_med <- gl.drop.pop(bass_all_incl_med, pop.list=c("sp_Catalonia","sp_Portugal_N","sp_Portugal_S"))

#check pops again - should have removed three pops from above
levels(bass_all_excl_med$pop)

#use dartR to make fst matrix. 1000 bootstraps
bass_all_excl_med_fst<-gl.fst.pop(bass_all_excl_med, nboots=1000)

#extract pairwise fsts
fst_pairwise_all_pops<-bass_all_excl_med_fst$Fsts
#extract p values
pval_pairwise_all_pops<-bass_all_excl_med_fst$Pvalues

#get pops
areas_keep_all<-as.vector(unique(bass_all_excl_med$pop))

#give matrix size
n=35
# Initialize a matrix to store the reordered lower triangular part with NAs
reordered_matrix_all <- matrix(NA, nrow = n, ncol = n, dimnames = list(areas_keep_all, areas_keep_all))


# Reorder the lower triangular part while preserving NAs
for (i in 1:n) {
  for (j in 1:i) {
    if (!is.na(fst_pairwise_all_pops[areas_keep_all[i], areas_keep_all[j]])) {
      reordered_matrix_all[i, j] <-fst_pairwise_all_pops[areas_keep_all[i], areas_keep_all[j]]
    } else {
      reordered_matrix_all[i, j] <- fst_pairwise_all_pops[areas_keep_all[j], areas_keep_all[i]]
    }
  }
}

#convert to data frame

#convert fsts to data frame
fst.mat_test_all <- reshape::melt(reordered_matrix_all) %>%
  mutate(value = as.numeric(value))

###do p values reorganisation

#give matrix size
n=35
# Initialize a matrix to store the reordered lower triangular part with NAs
reordered_matrix_all_pval <- matrix(NA, nrow = n, ncol = n, dimnames = list(areas_keep_all, areas_keep_all))


# Reorder the lower triangular part while preserving NAs for pvals
for (i in 1:n) {
  for (j in 1:i) {
    if (!is.na(pval_pairwise_all_pops[areas_keep_all[i], areas_keep_all[j]])) {
      reordered_matrix_all_pval[i, j] <-pval_pairwise_all_pops[areas_keep_all[i], areas_keep_all[j]]
    } else {
      reordered_matrix_all_pval[i, j] <-pval_pairwise_all_pops[areas_keep_all[j], areas_keep_all[i]]
    }
  }
}


#convert pvalues to data frame

fst.mat_test_pvals <- reshape::melt(reordered_matrix_all_pval) %>%
  mutate(value = as.numeric(value))

#get no. pvals in data set
tests_all<-(nrow(fst.mat_test_pvals)-sum(is.na(fst.mat_test_pvals$value)))


# Adjust p-values and create significance labels in one step
fst.mat_test_all_pvals <- fst.mat_test_pvals %>%
  mutate(
    adj = p.adjust(value, method = "BH", n = tests_all),
    stars = cut(adj, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf), 
                labels = c("a", "b", "c", ""))
  )






# Convert minus values to zero
fst.mat_test_all$value[fst.mat_test_all$value < 0] = 0


# Fst italic label
fst.label = expression(italic("F")[ST])

#plot heatmap all seasons, exl med ices rectangles.

fst_all<-ggplot(data = fst.mat_test_all, aes(x = X2, y = X1, fill = value))+
  geom_tile(color = "black")+ scale_fill_gradient(low = "white", high = "red", name=fst.label, na.value = "white") + 
  geom_text(label=fst.mat_test_all_pvals$stars,color="black", size=2) +
  #ggtitle(expression(atop("Pairwise FST all spawning and feeding samples, WC (1984) Dartr")))+
  ggtitle(bquote(Pairwise~F[ST]~spawning~and~feeding~samples~from~around~UK ))+
  labs( x = "ICES rectangle", y = "ICES rectangle") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 8, hjust = 1),axis.text.y = element_text(size = 8)) + 
  coord_fixed()


pdf(file="output/fst_analysis/fst_all_excl_med_sep_24.pdf")
fst_all
dev.off()
