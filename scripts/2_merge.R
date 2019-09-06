#merge files e.g. folder tile10
#first ensure you are looking to the rigth data, list of .tif files

#running from the terminal this line below, you should list all the tif files in each folder that 
#called with a name starting bfast_ and in particular each tif file in these folder that end as _2019.tif
#(to locate the raw data instead the trheshold data)

#echo `ls bfast_*/*_2019.tif`; echo

#if the list I get from the command line is the correct, I could run the following

#setworking directory to the directory I am going to work on the merging the files 
getwd()

setwd("~/Sudan_BFAST/data/laura_bfast_sudan/")
#list files in the folder
list.files(path = ".", pattern = NULL, all.files = FALSE,
           full.names = FALSE, recursive = FALSE,
           ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)

#list the files I want to merge (raw data)
system("
       for r in `ls bfast_*/*_2019.tif`; do
       echo $r
       done
       ")

system("gdal_merge.py -v -co \"COMPRESS=LZW\" -co \"BIGTIFF=YES\" -o \"tmp_tile10_merge.tif\"  $(ls bfast_*/*_2019.tif)")
#check again the list of tiles, you will see a new file called "tmp_tile10_merge.tif"
#now you can eventually repeat the same steps for the threshold tif files

#list the files I want to merge (threshold.tif data)
system("
       for r in `ls bfast_*/*_threshold.tif`; do
       echo $r
       done
       ")
system("gdal_merge.py -v -co \"COMPRESS=LZW\" -co \"BIGTIFF=YES\" -o \"tmp_tile10_threshold_merge.tif\"  $(ls bfast_*/*_threshold.tif)")


