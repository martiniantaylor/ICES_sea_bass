# ICES_sea_bass
R scripts supporting the analyses in Taylor et al. (2025), ICES Journal of Marine Science, on the population genetics of European sea bass (Dicentrarchus labrax).

## Script: 04_basic_stats_v2.R

### Overview 
This script calculates basic genetic diversity statistics - observed heterozygosity (Ho), expected heterozygosity (He), and inbreeding coefficient (Fis) for spawning and feeding populations of sea bass. The results are formatted into publication-ready tables using gt.

### Required Files

- **Genlight files:**
The script requires the following genotype datasets: 

- `data/spawning_incl_med.gl` (41022 loci, 343 individuals, genlight format)
- `data/bass_reduced_ices_feeding.gl` (41022 loci, 558 individuals, genlight format)
- **ICES area data:**
- `data/ices_area_rectangle.csv` (Used to group samples in ices rectangle format by ICES region)


### Output
The script generates tables summarizing genetic diversity statistics for spawning and feeding populations. These tables are formatted for inclusion in the manuscript's supplementary materials.

***

## Script: 05_pairwise_fst_analysis_v4.R

### Overview 
This script performs pairwise F<sub>ST</sub> analysis on genetic data from bass populations across ICES rectangles and ICES areas, including independent analysis of spawning and feeding seasons.


### Required Files

- **Genlight files:**
The script requires the following genotype datasets: 
-  `data/spawning_incl_med.gl` (41022 loci, 343 individuals, genlight format)
-  `data/bass_reduced_ices_feeding.gl` (41022 loci, 558 individuals, genlight format))

### Output

Heatmaps of pairwise F<sub>ST</sub> for spawning and feeding samples and across all samples
Boxplot comparing spawning vs feeding F<sub>ST</sub>

***

## Script: 05c_overall_fst_v2.R

### Overview  
This script calculates bootstrapped overall *F*<sub>ST</sub> values for European sea bass (*Dicentrarchus labrax*) using the `finepop2` package. The script processes genetic data for spawning and feeding populations in both the Atlantic and Mediterranean. It converts genetic data into the `genepop` format, extracts population names, and computes global *F*<sub>ST</sub> estimates along with standard errors.  

### Required Files 
- **Genlight files:**
The script requires the following genotype datasets:  
- `data/bass_reduced_ices_spawning.gl` (spawning populations dataset in genlight format)  
- `data/bass_reduced_ices_feeding.gl` (feeding populations dataset in genlight format)  
- `data/spawning_incl_med.gl` (spawning populations including Mediterranean samples)  
- **ICES area data:**
- `data/ices_area_rectangle.csv` (Used to group samples in ices rectangle format by ICES region)
Additionally, it requires the following external scripts:  
- `scripts/jenkins_genepop_func.R` (R function to export a genind object in genepop format)  from [Tom Jenkins' repo](https://github.com/Tom-Jenkins/utility_scripts/tree/master)

### Output
Bootstrapped overall bootstrapped overall *F*<sub>ST</sub> values for various geographic sample groupings.


***

## Script: 05d_pop_specific_fst_v2.R

### Overview 
This script calculates **population-specific FST** values using **FinePop2** and generates a combined plot for **feeding and spawning** data, including populations from **Portugal and the Mediterranean**. The results contribute to **Figure 4** in the sea bass manuscript.

### Required Files
The script requires the following input files:

- **Genlight files:**
  - `data/spawning_incl_med.gl` (Spawning + Mediterranean populations, 41,022 loci, 343 individuals, genlight format))
  - `data/bass_reduced_ices_feeding.gl` (Feeding dataset, 41,022 loci, 558 individuals, genlight format))
- **ICES area data:**
  - `data/ices_area_rectangle.csv` (Used to group samples in ices rectangle format by ICES region)

Additionally, the script sources a  helper function:

- `jenkins_genepop_func.R`: (R function to export a genind object in genepop format) from [Tom Jenkins' repo](https://github.com/Tom-Jenkins/utility_scripts/tree/master)

### Output
Figure 4 in manuscript. Combined plot of feeding and spawning population specific FST

***

## Script: 06_forceplots_fst_v3.R

### Overview 
This script generates forceplots of pairwise Fst data for feeding and spawning samples of sea bass. The script reads in pairwise Fst matrices, processes the data, and outputs the figures to PDF files.

### Required Files

The following files are required to run the script:

1. `data/reordered_matrix_spawn.txt`: Pairwise Fst matrix for spawning samples.
2. `data/reordered_matrix_feed_fst.txt`: Pairwise Fst matrix for feeding samples.
3. `data/ices_area_rectangle.csv`: CSV file containing the ices area for each ices rectangle.

### Output

The script generates the following output files:

1. `output/fst_analysis/test_force_sp.pdf`: Forceplot for spawning samples.
2. `output/fst_analysis/test_force_fe.pdf`: Forceplot for feeding samples.

***

## Script: 08_mantel_tests_v2.R 

### Overview 
This script performs Mantel tests to assess the correlation between pairwise genetic differentiation (Fst) and pairwise geographic distance for feeding and spawning populations. The results are visualized using scatter plots with linear regression lines.

### Required Files
The following files must be present for the script to run correctly:

1. **Fst Matrices** (generated and saved in the 05_pairwise_fst script):
   - `data/reordered_matrix_spawn.txt`
   - `data/reordered_matrix_feed_fst.txt`

2. **Geographic Distance Matrices** (generated and saved in the pairwise geographic distance script):
   - `output/sea_distances/Contour_PairwiseMarineDistances_feeding.txt`
   - `output/sea_distances/Contour_Pairwise_Marine_Distances_spawning.txt`

### Output
- Mantel test results printed to the console.
- Two PDF plots showing the correlation between genetic distance and geographic distance.

***

## Script: 09_pairwise_sea_distance_v2.R 

### Overview 
This script calculates pairwise sea distances between sea bass spawning and feeding sites using a marine distance function. The script relies on shapefiles and a function sourced from Jorge Assis's GitHub repository (https://raw.githubusercontent.com/jorgeassis/marineDistances/master/Script.R).

### Required Files

The following files are required to run the script: shapefiles from [NGDC GSHHG dataset](https://www.ngdc.noaa.gov/mgg/shorelines/data/gshhg/latest/) 

1. `data/map_shape_files/h/GSHHS_h_L1.shp`: Shapefile for global coastlines.
2. `data/map_shape_files/h/GSHHS_h_L1.dbf`: Database file for the shapefile.
3. `data/map_shape_files/h/GSHHS_h_L1.shx`: Index file for the shapefile.
4. `data/sea_distance/spawn_sites.csv`: CSV file containing latitude and longitude data for spawning sites.
5. `data/sea_distance/feed_sites.csv`: CSV file containing latitude and longitude data for feeding sites.

### Output

The script generates the following output files:

1. `output/sea_distances/Contour_spawning_sites.txt`: Pairwise sea distances for spawning sites.
2. `output/sea_distances/Contour_feeding_sites.txt`: Pairwise sea distances for feeding sites.

***

## Script: 10_pca_v3.R 

### Overview 
This script performs Principal Component Analysis (PCA) on sea bass spawning data, including Atlantic, Mediterranean, and Portugal populations. The script generates PCA plots and computes the variance explained by the principal components for Figure 5 of manuscript.

### Required Files

The following files are required to run the script:

1. `data/spawning_incl_med.gl`: Genlight file containing the genetic data for spawning populations including Atlantic, Mediterranean, and Portugal.
2. `data/bass_reduced_ices_spawning.gl`: Genlight file containing the genetic data for Atlantic spawning populations.
3. `data/ices_area_rectangle.csv`: CSV file containing the ICES area for each ICES rectangle.

### Output

The script generates the following output files:

1. `output/pca/pca_med_cat.pdf`: PCA plot including Mediterranean populations (Figure 5a).
2. `output/pca/pca_spawning_atl_ices_area.pdf`: PCA plot for Atlantic populations (Figure 5b).

***

## Script:  11_admixture_distance_plot_V3.R 

### Overview 
This script generates plots of admixture proportions against geographic distance from the Mediterranean for sea bass. The script reads in precomputed pairwise sea distances and genetic data, processes the data, and outputs the figures to PDF files.

### Required Files

The following files are required to run the script:

1. `output/sea_distances/Contour_spawning_incl_med_18_pops.txt`: Pairwise sea distances matrix.
2. `data/spawning_incl_med.gl`: Genetic data for spawning and Mediterranean samples, genlight format)
3. `data/ices_area_rectangle.csv`: CSV file containing the ICES area for each ICES rectangle.

### Output

The script generates the following output files:

1. `output/admixture/spawning_dist_zoom_sept24.pdf`: Plot of admixture proportions against distance.

***

## Script: 12a_pca_adapt_v1.R 

### Overview 
This script identifies outliers in spawning populations using PCAadapt and saves the identified outliers as a new genlight file. The script processes genetic data, performs PCA-based adaptation analysis, and identifies statistically significant outliers.

### Required Files

The following files are required to run the script:

1. `data/outflank/spawning_incl_med.gl`: genlight file containing the genetic data for spawning populations including Mediterranean samples.

### Output

The script generates the following output files:

1. `output/pca_adapt/sp_outlying_loci_pca_adapt.gl`: Genlight file containing the identified outlier SNPs.
2. `pcaasapt_outliers_vec`: File containing the names of the outlier loci.

***

## Script:  12b_bayescan_V2.R 

### Overview 

This script uses BayeScan to detect outliers in SNP data from Atlantic and Mediterranean sea bass populations. The script processes the data, runs BayeScan, performs convergence diagnostics, and visualizes the results.

### Required Files

The following files are required to run the script:

1. `data/outflank/spawning_incl_glport.gl`: genlight file containing the genetic data for Atlantic spawning populations including Portugal N.
2. `data/outflank/spawning_incl_med.gl`: genlight file containing the genetic data for Mediterranean populations.

### Output

The script generates the following output files:

1. `output/bayescan/bayescan_sp_atl.txt`: BayeScan input file for Atlantic populations.
2. `output/bayescan/bayescan_sp_med.txt`: BayeScan input file for Mediterranean populations.
3. `output/bayescan/outlying_loci_bayescan.gl`: Genlight file containing the identified outlier SNPs.

***

## Script: 12c_outflank_v1.R 

### Overview 

This script uses OutFLANK to identify outliers in SNP data from Atlantic and Mediterranean sea bass populations. The script processes the data, runs OutFLANK, and visualizes the results.

### Required Files

The following files are required to run the script:

1. `data/spawning_incl_glport.gl`: genlight file containing the genetic data for Atlantic spawning populations (excluding Portugal S).
2. `data/spawning_incl_med.gl`: genlight file containing the genetic data for Mediterranean populations.

### Output

The script generates the following output files:

1. `output/outflank/he_fst_med.pdf`: Plot of He versus Fst for Mediterranean populations.
2. `output/outflank/he_fst_atl.pdf`: Plot of He versus Fst for Atlantic populations.
3. `output/pca_adapt/sp_outlying_loci_outflank.gl`: Genlight file containing the identified outlier SNPs.

***

## Script: 12d_combine_outlier_methods_V2.R 

### Overview 

This script combines SNPs identified as outliers from three different selection detection methods (OutFLANK, BayeScan, and PCAadapt). The script merges the outliers into a single list and identifies SNPs common across all methods.

### Required Files

The following files are required to run the script:

1. `output/bayescan/outlying_loci_bayescan.gl`: Genlight file containing outliers detected by BayeScan.
2. `output/outflank/sp_outlying_loci_outflank.gl`: Genlight file containing outliers detected by OutFLANK.
3. `output/pca_adapt/sp_outlying_loci_pca_adapt.gl`: Genlight file containing outliers detected by PCAadapt.
4. `data/outflank/spawning_incl_med.gl`: Genlight file containing the full genetic dataset for spawning populations including Mediterranean samples.

### Output

The script generates the following output files:

1. `output/outflank/sp_outlying_loci_all.gl`: Genlight file containing all unique outlier SNPs identified by different methods.
2. `output/outflank/sp_common_outliers_all_methods.txt`: Text file containing the names of SNPs identified as outliers by all methods.

***

## Script: 13_admixture_outliers_distance_v1.R 

### Overview 

This script generates admixture proportion plots against geographic distance from the Mediterranean for SNPs identified as outliers by three different selection detection methods (OutFLANK, BayeScan, and PCAadapt). The script processes the data, runs LEA for Mediterranean ancestry, and plots the results.

### Required Files

The following files are required to run the script:

1. `output/bayescan/outlying_loci_bayescan.gl`: Genlight file containing outliers detected by BayeScan.
2. `output/outflank/sp_outlying_loci_outflank.gl`: Genlight file containing outliers detected by OutFLANK.
3. `output/pca_adapt/sp_outlying_loci_pca_adapt.gl`: Genlight file containing outliers detected by PCAadapt.
4. `output/outflank/sp_outlying_loci_all.gl`: Genlight file containing all unique outlier SNPs identified by different methods.
5. `output/sea_distances/Contour_spawning_incl_med_18_pops.txt`: Pairwise sea distances matrix.
6. `data/ices_area_rectangle.csv`: CSV file containing the ICES area for each ICES rectangle.

### Output

The script generates the following output files:

1. `data/LEA_input/outlier_snps_all_df.geno`: GENO file for all outliers - this is the required format for LEA
2. `data/LEA_input/bayescan_snps_df.geno`: GENO file for BayeScan outliers - this is the required format for LEA
3. `data/LEA_input/outflank_snps_df.geno`: GENO file for OutFLANK outliers - this is the required format for LEA
4. `data/LEA_input/pcaadapt_snps_df.geno`: GENO file for PCAadapt outliers - this is the required format for LEA
5. `output/admixture/outliers_dist_all.pdf`: Plot of admixture proportions against distance for all outliers.
6. `output/admixture/outliers_dist_bayes.pdf`: Plot of admixture proportions against distance for BayeScan outliers.
7. `output/admixture/outliers_dist_outflank.pdf`: Plot of admixture proportions against distance for OutFLANK outliers.
8. `output/admixture/outliers_dist_pcaadapt.pdf`: Plot of admixture proportions against distance for PCAadapt outliers.
9. `output/admixture/model_stats.txt`: Text file containing the model statistics for each set of outliers.

***

## Script:  14_mapping_script_V4.R 

### Overview 

This script generates maps of sea bass feeding and spawning samples, including an inset map of southern samples. The script loads latitude and longitude data for the samples, as well as ICES rectangles, and produces maps with the specified data.

### Required Files

The following files are required to run the script:

1. `data/metadata/bass_samples_feeding_10.csv`: Latitude and longitude data for feeding samples.
2. `data/metadata/bass_samples_spawning_10.csv`: Latitude and longitude data for spawning samples.
3. `data/metadata/bass_samples_southern.csv`: Latitude and longitude data for southern samples.
4. `data/map_shape_files/ICES_Areas_20160601_cut_dense_3857.shp`: Shapefile for ICES areas.
5. `data/map_shape_files/ICES_Areas_20160601_cut_dense_3857.dbf`: Database file for ICES areas shapefile.
6. `data/map_shape_files/ICES_Areas_20160601_cut_dense_3857.shx`: Index file for ICES areas shapefile.

### Output

The script generates the following output files:

1. `output/maps/spawning_map.pdf`: Map of spawning samples.
2. `output/maps/feeding_map.pdf`: Map of feeding samples.
3. `output/maps/spawning_map_southern.pdf`: Combined map of spawning samples with southern inset.



