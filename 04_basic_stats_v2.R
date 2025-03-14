######################################################
#Script to output basic stats Ho, He etc to make the 
#Supp tables 2a&b in sea bass manuscript
#
#Separate tables for feeding and spawning samples
#
#
#Martin Taylor (martin.taylor@uea.ac.uk)
######################################################

# Load necessary libraries
library(dartR) #basic pop gen stats
library(dplyr) #data wrangling
library(gt) #output and manipulation in table format

#import saved spawning + MED file  - 41022 loci 343 individuals
bass_reduced_ices_area_spawning_MED<-gl.load("data/spawning_incl_med.gl")

#import feeding data set gl - 41022 loci, 558 indivs
bass_reduced_ices_feeding<-gl.load("data/bass_reduced_ices_feeding.gl")

#check pop levels
levels(bass_reduced_ices_area_spawning_MED$pop) #15 northeast atlantic pops + portugal (x2) and med
levels(bass_reduced_ices_feeding$pop) #20 populations

##get stats using dartR function for spawning pops
he_spawn<-gl.report.heterozygosity(
  bass_reduced_ices_area_spawning_MED,
  method = "pop",
  n.invariant = 0,
  plot.out = F,
  plot_theme = theme_dartR(),
  plot_colors_pop = discrete_palette,
  plot_colors_ind = two_colors,
  save2tmp = FALSE,
  verbose = NULL
)

#calculate percentage of polymorphic loci
he_spawn$poly_perc<-(1-he_spawn$monoLoc/he_spawn$polyLoc)*100
he_spawn$nInd<-round(he_spawn$nInd,0)

#select relevant columns from output 
spawn_he_table<-he_spawn%>% dplyr::select(1,2,8,12,18,19)
colnames(spawn_he_table)<-c("Ices rectangle", "N", "Ho", "He", "Fis", "% loci polymorphic")

#make outputs 3 decimal places
spawn_he_table_dp3<-format(spawn_he_table, digits = 3, row.names=FALSE)

# using gt to make tables....
spawn_stats<-spawn_he_table_dp3 |>
  gt() |>
  tab_header(
    title = "Bacic statistics for spawning populations"
  )|>
  fmt_number(
    columns =c(3:5),
    rows = everything(),
    decimals = 3,
  )|>
  fmt_number(
    columns =c(6),
    rows = everything(),
    decimals = 1,
  )

#save tables to pdf and word
spawn_stats |> gtsave("output/basic_stats/spawn_he_gttable_2025.pdf")
spawn_stats |> gtsave("output/basic_stats/spawn_he_gttable_2025.docx")


#get stats using dartR function for feeding samples
he_feed<-gl.report.heterozygosity(
  bass_reduced_ices_feeding,
  method = "pop",
  n.invariant = 0,
  plot.out = F,
  plot_theme = theme_dartR(),
  plot_colors_pop = discrete_palette,
  plot_colors_ind = two_colors,
  save2tmp = FALSE,
  verbose = NULL
)
#calculate % polymorphic loci
he_feed$poly_perc<-(1-he_feed$monoLoc/he_feed$polyLoc)*100
he_feed$nInd<-round(he_feed$nInd,0)

#get relevant output columns
feed_he_table<-he_feed%>% dplyr::select(1,2,8,12,18,19)
colnames(feed_he_table)<-c("Ices rectangle", "N", "Ho", "He", "Fis", "% loci polymorphic")

#make outputs 3 decimal places
feed_he_table_dp3<-format(feed_he_table, digits = 3, row.names=FALSE)


#using gt to make tables....
feed_stats<-feed_he_table_dp3 |>
  gt() |>
  tab_header(
    title = "Bacic statistics for feeding populations"
  )|>
  fmt_number(
    columns =c(3:5),
    rows = everything(),
    decimals = 3,
  )|>
  fmt_number(
    columns =c(6),
    rows = everything(),
    decimals = 1,
  )

#save tables to pdf and word
feed_stats |> gtsave("output/basic_stats/feed_he_gttable_2025.pdf")
feed_stats |> gtsave("output/basic_stats/feed_he_gttable_2025.docx")
