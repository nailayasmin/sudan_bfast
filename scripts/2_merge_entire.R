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

#the files from the list I need to pick up are (final_results_etc...) and when merging the threshold data, the file ending with
#"_threshold.tif"
#same for tile 00, 02,03, 04, 05, 06, 07

#for folder 08 the structure of the subfolder is:

# setwd("~/Sudan_BFAST/data/laura_bfast_sudan/")
#  list.files(path = ".", pattern = NULL, all.files = FALSE,
#             full.names = FALSE, recursive = FALSE,
#             ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)

#echo `ls bfast_*/*_2019.tif` > in the command line

#if the list I get from the command line is the correct, I could run the following
#list the files I want to merge (raw data)
# system("
#        for r in `ls *.tif`; do
#        echo $r
#        done
#        ")
# 
# system("gdal_merge.py -v -co \"COMPRESS=LZW\" -co \"BIGTIFF=YES\" -o \"final_results_tile08_merge_rawnaila.tif\"  $(ls *.tif)")
#now you can eventually repeat the same steps for the threshold tif files

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
getwd()
#https://unix.stackexchange.com/questions/60849/find-files-in-multiple-folder-names
#https://unix.stackexchange.com/questions/60849/find-files-in-multiple-folder-names/390333#390333

system("
       for r in `find . -type f -name final_results*.tif `; do
       echo $r
       done
       ")
#find . -type f -name "final_results*.tif"

system("gdal_merge.py -v -co \"COMPRESS=LZW\" -co \"BIGTIFF=YES\" -o \"raw_final_results_alltiles.tif\"  $(find . -type f -name final_results*.tif)")

