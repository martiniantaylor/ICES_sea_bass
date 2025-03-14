
######################################################
#Script to calculate pairwise sea distances 
#needs file of lat lons and pops
#relies on marine distance function by jorgeassis
#
#note-the shp file needs to be with the dbf and shx files
#download the shape files from here 
#(https://www.ngdc.noaa.gov/mgg/shorelines/data/gshhg/latest/)
#
#
#Martin Taylor (martin.taylor@uea.ac.uk)
######################################################


#source the marine distnaces function
source("https://raw.githubusercontent.com/jorgeassis/marineDistances/master/Script.R")

#load shapefiles - downloaded above
global.polygon <- "data/map_shape_files/h/GSHHS_h_L1.shp"

## Run the function - saves outputs to file. Lons and lats saved in data folder
contour( global.polygon = global.polygon , file= "data/sea_distance/spawn_sites.csv" , file.sep = "," , file.dec = "." , file.strucutre = 1 , file.header = TRUE , resolution = 0.01 , buffer = c(1,1,1,1) , export.file = TRUE )
contour( global.polygon = global.polygon , file= "data/sea_distance/feed_sites.csv" , file.sep = "," , file.dec = "." , file.strucutre = 1 , file.header = TRUE , resolution = 0.01 , buffer = c(1,1,1,1) , export.file = TRUE )

## file : the main file with the locations; should be text delimited
## global.polygon: the path of the polygon
## file.strucutre: the main file structure: 1 to “Name Lon Lat” or 2 to “Name Lat Lon”
## file.header: define if the text file has a header with the column names (TRUE or FALSE)
## resolution: the resolution of the study area and the buffer to use around the sites. 
## buffer: the buffer can be a simple value or a vector such as c(xmin,xmax,ymin,ymax). 
## export.file: file to export the results as a text delimited file (TRUE or FALSE)

