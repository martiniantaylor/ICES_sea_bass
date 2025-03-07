######################################################
#Population specific fst script using Finepop2
#Figure 4 in sea bass manuscript
#This will make combined plot using both feeding 
#and spawning data incl Portugal and Med pops

#Martin Taylor (martin.taylor@uea.ac.uk)
######################################################

#load libs
library(FinePop2)
library(dartR)

#these are required scripts
source("scripts/jenkins_genepop_func.R")
#This is modified to make it work with generated genepop files. The way markers are counted is modified.
source("scripts/finepoptest.R")

#import saved spawning + MED file  - 41022 loci 343 individuals
bass_reduced_ices_area_spawning_MED<-gl.load("data/spawning_incl_med.gl")

#import feeding data set gl - 41022 loci, 558 indivs
bass_reduced_ices_feeding<-gl.load("data/bass_reduced_ices_feeding")

#look at pop order in genlight
levels(bass_reduced_ices_area_spawning_MED$pop)
levels(bass_reduced_ices_feeding$pop)

#sort genlight objects by population
bass_reduced_ices_area_spawning_MED_sort<-gl.sort(bass_reduced_ices_area_spawning_MED, sort.by = "pop", verbose = NULL)
bass_reduced_ices_feeding_sort<-gl.sort(bass_reduced_ices_feeding_sort, sort.by = "pop", verbose = NULL)

#rename pops as _sp and _fe to allow separation later on into feeding and spawning
names_spawn <- data.frame(pop = unique(bass_reduced_ices_area_spawning_MED_sort$pop)) %>%
  mutate(new = paste("sp", pop, sep = "_"))

names_feed <- data.frame(pop = unique(bass_reduced_ices_feeding_sort$pop)) %>%
  mutate(new = paste("fe", pop, sep = "_"))

#remove column names
colnames(names_spawn)=NULL
colnames(names_feed)=NULL

#save as csv files
write.table(names_spawn, file="data/temp_files/recode_pops_spawn_2025.csv", row.names = F, col.names = F, sep =",")
write.table(names_feed, file="data/temp_files/recode_pops_feed_2025.csv", row.names = F, col.names = F, sep =",")

#rename populations using dartr function
bass_reduced_ices_area_spawning_MED_sort_rename<-gl.recode.pop(bass_reduced_ices_area_spawning_MED_sort, pop.recode="data/temp_files/recode_pops_spawn_2025.csv", recalc = TRUE, mono.rm = FALSE, verbose = 5)
bass_reduced_ices_feeding_sort_rename<-gl.recode.pop(bass_reduced_ices_feeding_sort, pop.recode="data/temp_files/recode_pops_feed_2025.csv", recalc = TRUE, mono.rm = FALSE, verbose = 5)

#check pop levels
levels(bass_reduced_ices_area_spawning_MED_sort_rename$pop)
levels(bass_reduced_ices_feeding_sort_rename$pop)

#combine feeding and spawning
combined_feed_sp<-rbind(bass_reduced_ices_area_spawning_MED_sort_rename, bass_reduced_ices_feeding_sort_rename)

#save file for other scripts
gl.save(combined_feed_sp, file="data/combined_incl_med.gl", verbose=NULL)

#check pops
levels(combined_feed_sp$pop)

#make genind object from genlight
combined_feed_sp_genind<-gl2gi(combined_feed_sp)

#check pops are still as expected
levels(combined_feed_sp_genind$pop)

#convert to genepop
genind2genepop_jenk(combined_feed_sp_genind, file="data/finepop2/combined_feed_sp_2025.txt")

#put pop names into vector
pop_names_all<-as.character(unique(combined_feed_sp_genind$pop))

#save pop names to file as vector
sink("data/finepop2/pop_names_all.txt")
cat(pop_names_all)
sink()

#import genepop file and populations file for finepop
combined_finepop<-read.GENEPOP(genepop="data/finepop2/combined_feed_sp_2025.txt", popname="data/finepop2/pop_names_all.txt")

#check popnames (again!)
combined_finepop$pop_names

#FST estimation of feeding and spawning incl med----------------
combined_finepop.popspFST<-pop_specificFST(combined_finepop) 
#print fst
print(combined_finepop.popspFST$fst)

#add names to object
combined_finepop.popspFST$names<-row.names(combined_finepop.popspFST$fst)


#organise outputs 

# Create a new data frame necessary data from fst outputs
data <- data.frame(
  fst = combined_finepop.popspFST$fst$FST,
  se = combined_finepop.popspFST$fst$SE,
  pop = row.names(combined_finepop.popspFST$fst)
)

# Extract color codes (assuming lowercase letters of length 1-2 represent them)
data <- data %>%
  mutate(colours = str_extract(pop, "[a-z]{1,2}"))

#this function removes the first part and separator of y-axis pop names
myfun = function(x) word(x, 2, sep = "_")

#make pure pop column in data using existing function above
data$ices_rect<-myfun(data$pop)

#rename sp and fe to Spawning and Feeding
data$facet <- ifelse(data$colours == "sp", "Spawning", "Feeding")

#get ices area for each ices rectangle (prev saved in script 03)
ices_tab_uniq<-read.csv("data/ices_area_rectangle.csv", sep=",")

#add ices areas to data
data2<-left_join(data, ices_tab_uniq)

# Create a combined column with ICES rectangle and area
data2 <- data2 %>%
  mutate(combined = paste0(ices_rect, " (", ices_area, ")"))

# Rename specific populations to remove prefixes
rename_map <- c("sp_Catalonia" = "Catalonia",
                "sp_Portugal_N" = "Portugal_N",
                "sp_Portugal_S" = "Portugal_S")

data2$combined <- ifelse(data2$pop %in% names(rename_map), rename_map[data2$pop], data2$combined)

# Assign facet categories based on 'colours'
data2 <- data2 %>%
  mutate(facet = case_when(
    colours == "sp" ~ "Spawning",
    TRUE ~ "Feeding"
  ))

# Ensure specific populations are set to "Ancestral"
ancestral_pops <- c("Catalonia", "Portugal_N", "Portugal_S")
data2$facet[data2$combined %in% ancestral_pops] <- "Ancestral"

#plot with ices areas added to labels
p3<-data2%>% 
  mutate(group = tidytext::reorder_within(combined, -fst, facet), fst) %>% 
  ggplot(aes(fst, group, colour=facet)) +
  geom_point() +
  geom_errorbar(aes(xmin=fst-(2*se), xmax=fst+(2*se), width=.2)) +
  labs(x = "Population specific Fst", y = "ICES rectangle", color="Season", title = "")+
  tidytext::scale_y_reordered(labels=)+
  facet_grid(~data2$facet ~ .,scales = "free", space = "free")+
  theme(strip.text.y = element_text(angle = 0), legend.position="none")

pdf(file="output/finepop/pop_specific_fst_facet_mar_25.pdf")
p3
dev.off()


