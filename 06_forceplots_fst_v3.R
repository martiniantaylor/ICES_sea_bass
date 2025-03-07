#ready for pub
#script to make force network plots

# Load required libraries
library(qgraph)
library(dplyr)

#get  populatons ices rectangles as a tibble
#spawning data
spawn_pops_tib <- tibble(ices_rect = pops_keep_spawn)
#feeding data
feed_pops_tib  <- tibble(ices_rect = pops_keep_feed)

#get ices area for each ices rectangle (prev saved in script 03)
ices_tab_uniq<-read.csv("data/ices_area_rectangle.csv", sep=",")

#combine rectangle with ices areas and convert iceas area to factor to allow colors 

#spawning
groups_sp <- spawn_pops_tib %>%
  left_join(ices_tab_uniq, by = "ices_rect") %>%
  distinct() %>%
  mutate(ices_area = factor(ices_area))

#feeding
groups_f <- feed_pops_tib %>%
  left_join(ices_tab_uniq, by = "ices_rect") %>%
  distinct() %>%
  mutate(ices_area = factor(ices_area))

# Load pairwise Fst matrices (saved in pairwise fst script)
Fst_matrix_spawn<-as.dist(read.table(file="data/reordered_matrix_spawn.txt"))
Fst_matrix_feed<-as.dist(read.table(file="data/reordered_matrix_feed_fst.txt"))

#convert -ve fsts to zero
Fst_matrix_spawn[Fst_matrix_spawn < 0] = 0
Fst_matrix_feed[Fst_matrix_feed < 0] = 0

#find maximum fst across feeding and spawning
max_fst<-max(Fst_matrix_spawn, Fst_matrix_feed)
print(max_fst)

#normalise and convert to similarity matrices to allow 
#direct comparision of line widths 

sim_sp <- 1 - (Fst_matrix_spawn / max_fst)
sim_fe <- 1 - (Fst_matrix_feed / max_fst)

# Debugging: Check similarity range
cat("Similarity matrix range (Spawn):", min(sim_sp), "-", max(sim_sp), "\n")
cat("Similarity matrix range (Feed):", min(sim_fe), "-", max(sim_fe), "\n")


# Define colors for plots
MyPalette_F <- c(
  "3.a.20" = "#DB9D85", "4.b" = "#C2A968", "4.c" = "#9DB469", 
  "6.a" = "#6DBC86", "7.a" = "#3DBEAB", "7.b" = "#4CB9CC", 
  "7.d" = "#87AEDF", "77.f" = "#BB9FE0", "7.g" = "#DA95CC"
)

MyPalette_S <- c(
  "4.c" = "#9DB469", "7.a" = "#3DBEAB", "7.b" = "#4CB9CC", 
  "7.d" = "#87AEDF", "7.e" = "#E494AB", "7.f" = "#BB9FE0"
)


pdf("output/fst_analysis/force_graphs_for_fst_col_2025_circle.pdf")#, width=1000, height=800)
par(mfrow = c(1, 2))
qgraph(sim_sp, layout='spring', vsize=3, title = "FST spawning samples", labels=pops_keep_spawn,groups=groups_sp_fact, color=MyPalette_S, cut =0.95, esize=7, colFactor=10, filetype="pdf", filename="output/fst_analysis/test_force_sp.pdf")
qgraph(sim_fe, layout='spring', vsize=3, title = "Fst feeding samples", labels=pops_keep_feed, groups=groups_F_fact, color= MyPalette_F, cut=0.95,  esize=7, colFactor=10, filetype="pdf", filename="output/fst_analysis/test_force_fe.pdf")
dev.off()

