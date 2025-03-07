# Packages ----------------------------------------------------------------
#these are the libraries required for the sea bass analysis

#general
library(tidyverse)
library(mmod)
library(ggtext)#ggplot titles text manipulation

#stats scripts
library(pcadapt)
library(qvalue)

#data wrangling libs
library(reshape)
library(reshape2)
library(lubridate)
library(Matrix)
library(data.table)
library(gt)

#popgen libs
library(adegenet)
library(dartR)#general pop gen using genlight objects
library(LEA)#model based clustering
library(hierfstat)#overall fst

#gis ices type mapping libs
library(RstoxFDA)#ices rectangle and area info
library(mapplots)
library(lubridate)#messing with dates
library(rnaturalearth)#map data sets
library(sf)#plotting maps

#outliers
library(qvalue)
library(pcadapt)
library(OutFLANK)
