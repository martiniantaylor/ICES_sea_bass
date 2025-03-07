####get bootstrapped overall FST using finepop2 package

source("scripts/jenkins_genepop_func.R")
source("scripts/01_load_libraries_v1.R")
#this finepoptest script was modified by mit to change  "\b" into "\n" on line 16
#otherwise it doesn't work!
source("scripts/finepoptest.R")

#get data files
#spawning gl - 41022 loci, 286 indivs
bass_reduced_ices_spawning<-gl.load("data/bass_reduced_ices_spawning")

#import feeding data set gl - 41022 loci, 558 indivs
bass_reduced_ices_feeding<-gl.load("data/bass_reduced_ices_feeding")

#sort genlight object to get all pops together
bass_reduced_ices_spawning_sort<-gl.sort(bass_reduced_ices_spawning, sort.by = "pop", verbose = NULL)
bass_reduced_ices_feeding_sort<-gl.sort(bass_reduced_ices_feeding, sort.by = "pop", verbose = NULL)

#make genind object
bass_red_spawning.gi<-gl2gi(bass_reduced_ices_spawning_sort)
bass_red_feeding.gi<-gl2gi(bass_reduced_ices_feeding_sort)

#convert to genepop uses function from here: https://github.com/Tom-Jenkins/utility_scripts/blob/master/TJ_genind2genepop_function.R
#mit renamed the function in the script to genind2genepop_jenk

genind2genepop_jenk(bass_red_spawning.gi, file="data/bass_red_spawning.gi.txt")
genind2genepop_jenk(bass_red_feeding.gi, file="data/bass_red_feeding.gi.txt")

#get popnames
pop_names_sp<-as.character(unique(bass_red_spawning.gi$pop))
pop_names_fe<-as.character(unique(bass_red_feeding.gi$pop))

#save names to file as vector
sink("data/finepop2/spawn_atl_names.txt")
cat(pop_names_sp)
sink()

#save names to file as vector
sink("data/finepop2/feed_atl_names.txt")
cat(pop_names_fe)
sink()


#atlantic spawning - make sure to use the modified read.GENEPOP function from finepoptest.R above
popdata_sp<-read.GENEPOP(genepop="data/bass_red_spawning.gi.txt", popname="data/finepop2/spawn_atl_names.txt")

#atlantic feeding
popdata_fe<-read.GENEPOP(genepop="data/bass_red_feeding.gi.txt", popname="data/finepop2/feed_atl_names.txt")

#check pop sizes
popdata_sp$pop_sizes
popdata_fe$pop_sizes

#calc global fst and se
atl_spawn.globalFST<-globalFST(popdata_sp)
print(atl_spawn.globalFST)

#$fst
#[1] 0.0003161239
#
#$se
#[1] 5.463942e-05

#calc global fst and se
atl_feed.globalFST<-globalFST(popdata_fe_atl)
print(atl_feed.globalFST)

#$fst
#[1] 0.0002238124
#
#$se
#[1] 3.214723e-05


#get spawning pops + MED and Portugal samples

#import saved spawning + MED file 
bass_reduced_ices_area_spawning_MED<-gl.load("data/spawning_incl_med.gl")

#sort by population
bass_reduced_ices_area_spawning_MED_sort<-gl.sort(bass_reduced_ices_area_spawning_MED, sort.by = "pop", verbose = NULL)
#make genind
spawning_incl_med.sort.gi<-gl2gi(bass_reduced_ices_area_spawning_MED_sort)
#convert to genepop
genind2genepop_jenk(spawning_incl_med.sort.gi, file="data/bass_med_spawning.gi.txt")

#get pop names as vector
pop_names4<-as.character(unique(spawning_incl_med.sort.gi$pop))

#save names to file as vector
sink("data/finepop2/spawn_med_names.txt")
cat(pop_names4)
sink()

#med spawning
popdata_med<-read.GENEPOP(genepop="data/bass_med_spawning.gi.txt", popname="data/finepop2/spawn_atl_names.txt")

#calc global fst and se
med_spawn.globalFST<-globalFST(popdata_med)
print(med_spawn.globalFST)

#$fst
#[1] 0.006207072
#
#$se
#[1] 7.966732e-05


##########Global FST for ices areas#####################

#get ices area for each ices rectangle
#label by ices area-
ices_tab_uniq<-read.csv("data/ices_area_rectangle.csv", sep=",")

#get existing pops and join with ices area data
spawn_ices_rects <- tibble(ices_rect = bass_reduced_ices_spawning$pop) %>%
  left_join(ices_tab_uniq, by = "ices_rect")

#save as file with no colnames or row names
write.table(spawning_areas, file="data/temp_files/spawning_areas_rename.txt", row.names = F, col.names = F, quote = FALSE, sep =",")

#rename populations using dartr function
bass_reduced_ices_area_spawning<-gl.recode.pop(bass_reduced_ices_spawning, pop.recode="data/temp_files/spawning_areas_rename.txt", recalc = TRUE, mono.rm = FALSE, verbose = 5)
#sanity check 5 - make sure that this ices area and not rectangles still
bass_reduced_ices_area_spawning$pop

#sort genlight object
bass_reduced_ices_area_spawning_sort<-gl.sort(bass_reduced_ices_area_spawning, sort.by = "pop", verbose = NULL)

#make genind
bass_reduced_ices_area_spawning_sort.gi<-gl2gi(bass_reduced_ices_area_spawning_sort)
#convert to genepop
genind2genepop_jenk(bass_reduced_ices_area_spawning_sort.gi, file="data/bass_atl_icesarea_spawning.gi.txt")

#get pop names
pop_names5<-as.character(unique(bass_reduced_ices_area_spawning_sort.gi$pop))

#save names to file as vector
sink("data/finepop2/spawn_atl_area_names.txt")
cat(pop_names5)
sink()

#ices area spawning atlantic only - make sure finepoptest mod script is leaded see above
popdata_area_atl<-read.GENEPOP(genepop="data/bass_atl_icesarea_spawning.gi.txt", popname="data/finepop2/spawn_atl_area_names.txt")

#calc global fst and se
spawn.area_atl.globalFST<-globalFST(popdata_area_atl)
print(spawn.area_atl.globalFST)

#$fst
#[1] 0.0002235722
#
#$se
#[1] 4.142294e-05

#############atlantic feeding area global fst###################

#get existing pops and join with ices area data
feed_ices_rects <- tibble(ices_rect = bass_reduced_ices_feeding$pop) %>%
  left_join(ices_tab_uniq, by = "ices_rect")

#save as file with no colnames or row names
write.table(feeding_areas, file="data/temp_files/feeding_areas_rename.txt", row.names = F, col.names = F, quote = FALSE, sep =",")

#rename populations using dartr function
bass_reduced_ices_area_feeding<-gl.recode.pop(bass_reduced_ices_feeding, pop.recode="data/temp_files/feeding_areas_rename.txt", recalc = TRUE, mono.rm = FALSE, verbose = 5)
#sanity check - make sure that this iceas area and not rectangles still
bass_reduced_ices_area_feeding$pop

#sort genlight object
bass_reduced_ices_area_feeding_sort<-gl.sort(bass_reduced_ices_area_feeding, sort.by = "pop", verbose = NULL)

#make genind
bass_reduced_ices_area_feeding_sort.gi<-gl2gi(bass_reduced_ices_area_feeding_sort)
#convert to genepop
genind2genepop_jenk(bass_reduced_ices_area_feeding_sort.gi, file="data/bass_atl_icesarea_feeding.gi.txt")

#get pop names
pop_names6<-as.character(unique(bass_reduced_ices_area_feeding_sort.gi$pop))

#save names to file as vector
sink("data/finepop2/feed_atl_area_names.txt")
cat(pop_names6)
sink()

#ices area feeding atlantic only
popdata_area_atl_feed<-read.GENEPOP(genepop="data/bass_atl_icesarea_feeding.gi.txt", popname="data/finepop2/feed_atl_area_names.txt")

#calc global fst and se
feed.area_atl.globalFST<-globalFST(popdata_area_atl_feed)
print(feed.area_atl.globalFST)


#$fst
#[1] 0.0001030398
#
#$se
#[1] 0.00002616892