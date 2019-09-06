#merge files e.g. folder tile10 was an example in case you need to merge some files from a specific folder/ for a specific tile
#here you will repeat the same procedure to merge the raw data and the threshold data from all the tiles folders

#first of all, ensure the structure of the files and folder is homogeneous for each folder 
#if not, the list of files should be in line with the names and subfolders of the files I want to call/merge

#e.g. looking for the bfast_results01 folder:
#to get the files (tif files) I need to merge, I need to enter in the subfolder called

#ts_sd_2007_2019_tile_01_daietti_PARAM_O_1_H_ROC_T_OC_F_h_Overall_2007_2010_2019/

#for example
#setwd("~/4_sudan/data/bfast/tiles/bfast_results_sudan/bfast_results01/ts_sd_2007_2019_tile_01_daietti_PARAM_O_1_H_ROC_T_OC_F_h_Overall_2007_2010_2019/")
setwd("~/Sudan_BFAST/data/laura_bfast_sudan/")
list.files(path = ".", pattern = NULL, all.files = FALSE,
           full.names = FALSE, recursive = FALSE,
           ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)

#list the files I want to merge (threshold.tif data)
system("
       for r in `ls bfast_*/*_threshold.tif`; do
       echo $r
       done
       ")
system("gdal_merge.py -v -co \"COMPRESS=LZW\" -co \"BIGTIFF=YES\" -o \"th_final_results_tile08_merge.tif\"  $(ls bfast_*/*_threshold.tif)")

#####################################
###################################
#merge all the files: raw data are located in each bfast_result#tilefolder and subfolder that start with ts_sd

setwd("~/Sudan_BFAST/data/laura_bfast_sudan/")

system("
       for r in `find . -type f -name final_results*.tif `; do
       echo $r
       done
       ")
#find . -type f -name "final_results*.tif"

system("gdal_merge.py -v -co \"COMPRESS=LZW\" -co \"BIGTIFF=YES\" -o \"raw_final_results_alltiles.tif\"  $(find . -type f -name final_results*.tif)")

#read the merged products



sdn_bfast_threshold <- raster("~/Sudan_BFAST/data/bfast_result/th_final_results_alltiles_merge.tif")
sdn_bfast_mon_year <- raster("~/Sudan_BFAST/data/bfast_result/raw_final_results_alltiles.tif", band=1)
sdn_bfast_magnitude <- raster("~/Sudan_BFAST/data/bfast_result/raw_final_results_alltiles.tif",band=2)
sdn_bfast_error <- raster("~/Sudan_BFAST/data/bfast_result/raw_final_results_alltiles.tif",band=3)


# display results

plot(sdn_bfast_threshold)
plot(sdn_bfast_magnitude)
plot(sdn_bfast_mon_year)
plot(sdn_bfast_error)
plot(aoi,add=TRUE)

# Display defined values and years
plot(sdn_bfast_magnitude, zlim=c(-2000,-1))
plot(sdn_bfast_magnitude, zlim=c(0,0))
plot(sdn_bfast_mon_year, zlim=c(2015,2019))



