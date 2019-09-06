########################################
## Identifying patterns in BFAST output
########################################
## load parameters
#source('~/uga_activity_data/scripts/get_parameters.R')

# input files
#merged bfast name
bfastout <-paste0(bfast_dir,'rfl_raw_merged.tif')
index <- "savi"
# output file names
#bfast_mosaic <- "raw_merged_savi.tif"
#bfast_mosaic <- paste0(bfast_dir,  substr(basename(bfastout), 1, nchar(basename(bfastout))-4),'_savi.tif')
bfast_mosaic <- paste0(bfast_dir,  substr(basename(bfastout), 1, nchar(basename(bfastout))-4),'_', index,'.tif')
#result <- paste0(thres_dir,'2_raw_merged_savi_thf_mask.tif')
result <- bfast_mosaic 
forestmask <- paste0(lc_dir,'recl_all1234_1302frcl.tif')   #check raster(...)

## parameters
# factor to divide standard deviation
divide_sd <- 4

#You should have already the tif file with merged bfast results
### COMPRESS
system(sprintf("gdal_translate -co COMPRESS=LZW %s %s",
               paste0(bfast_dir_tiles,"raw_final_results_alltiles.tif"),
               bfast_mosaic
))
#############################################################
### CLEAN
#system(sprintf("rm %s",
#               paste0(bfast_dir_raw,"tmp_merge.tif")
#))
############################################

#check projections of bfast and forest mask
crs(raster(forestmask))
crs(raster(bfast_mosaic))
#################### reproject forest mask to latlong WGS84
fmask_wgs <- paste0(fo_dir,  substr(basename(forestmask), 1, nchar(basename(forestmask))-4),'_wgs.tif')
system(sprintf("gdalwarp -t_srs \"%s\" -overwrite -ot Byte -co COMPRESS=LZW %s %s",
               "EPSG:4326",
               forestmask,
             #paste0(fmask_dir,"tmp_proj.tif")
               fmask_wgs
))

crs(raster(fmask_wgs))
# clip fo mask to bfast output
fmask_wgs_clp <- paste0(fo_dir,  substr(basename(forestmask), 1, nchar(basename(forestmask))-4),'_wgs_clp.tif')
system(sprintf("gdal_translate -ot Byte -projwin %s %s %s %s -tr %s %s -co COMPRESS=LZW %s %s",
               extent(raster(bfast_mosaic))@xmin,
               extent(raster(bfast_mosaic))@ymax,
               extent(raster(bfast_mosaic))@xmax,
               extent(raster(bfast_mosaic))@ymin,
               res(raster(bfast_mosaic))[1],
               res(raster(bfast_mosaic))[2],
               
               #paste0(thres_dir,"tmp_proj.tif"),
               fmask_wgs,
               fmask_wgs_clp
))

# apply mask
system(sprintf("gdal_calc.py -A %s -B %s --B_band=2 --co COMPRESS=LZW --overwrite --outfile=%s --calc=\"%s\"",
               fmask_wgs_clp,
               bfast_mosaic,
               result,
               paste0("A*B")
))
# plot(raster(result))

## Post-processing ####
# calculate the mean, standard deviation, minimum and maximum of the magnitude band
# reclass the image into 10 classes
# 0 = no data
# 1 = no change (mean +/- 1 standard deviation)
# 2 = negative small magnitude change      (mean - 2 standard deviations)
# 3 = negative medium magnitude change     (mean - 3 standard deviations)
# 4 = negative large magnitude change      (mean - 4 standard deviations)
# 5 = negative very large magnitude change (mean - 4+ standard deviations)
# 6 = postive small magnitude change       (mean + 2 standard deviations)
# 7 = postive medium magnitude change      (mean + 3 standard deviations)
# 8 = postive large magnitude change       (mean + 4 standard deviations)
# 9 = postive very large magnitude change  (mean + 4+ standard deviations)
#################### SET NODATA TO NONE IN THE TIME SERIES STACK

tryCatch({
  
  outputfile   <- paste0(thres_dir,  substr(basename(bfast_mosaic), 1, nchar(basename(bfast_mosaic))-4),'_nomsk_threshold.tif')
  # r <- raster(result)
  # NAvalue(r) <- 0
  means_b2 <- cellStats( raster(result) , na.rm=TRUE, "mean") 
  mins_b2 <- cellStats( raster(result) , na.rm=TRUE,"min")
  maxs_b2 <- cellStats(  raster(result) ,na.rm=TRUE, "max")
  stdevs_b2 <- cellStats(  raster(result) ,na.rm=TRUE, "sd")/divide_sd
  system(sprintf("gdal_calc.py -A %s --co=COMPRESS=LZW --type=Byte --outfile=%s --calc='%s'
                 ",
                 result,
                 paste0(thres_dir,"tmp_threshold.tif"),
                 paste0('(A<=',(maxs_b2),")*",
                        '(A>',(means_b2+(stdevs_b2*4)),")*9+",
                        '(A<=',(means_b2+(stdevs_b2*4)),")*",
                        '(A>',(means_b2+(stdevs_b2*3)),")*8+",
                        '(A<=',(means_b2+(stdevs_b2*3)),")*",
                        '(A>', (means_b2+(stdevs_b2*2)),")*7+",
                        '(A<=',(means_b2+(stdevs_b2*2)),")*",
                        '(A>', (means_b2+(stdevs_b2)),")*6+",
                        '(A<=',(means_b2+(stdevs_b2)),")*",
                        '(A>', (means_b2-(stdevs_b2)),")*1+",
                        '(A>=',(mins_b2),")*",
                        '(A<', (means_b2-(stdevs_b2*4)),")*5+",
                        '(A>=',(means_b2-(stdevs_b2*4)),")*",
                        '(A<', (means_b2-(stdevs_b2*3)),")*4+",
                        '(A>=',(means_b2-(stdevs_b2*3)),")*",
                        '(A<', (means_b2-(stdevs_b2*2)),")*3+",
                        '(A>=',(means_b2-(stdevs_b2*2)),")*",
                        '(A<', (means_b2-(stdevs_b2)),")*2")
                 
  ))
  
}, error=function(e){})

####################  CREATE A PSEUDO COLOR TABLE
cols <- col2rgb(c("white","beige","yellow","orange","red","darkred","palegreen","green2","forestgreen",'darkgreen'))
pct <- data.frame(cbind(c(0:9),
                        cols[1,],
                        cols[2,],
                        cols[3,]
))

write.table(pct,paste0(thres_dir,"color_table.txt"),row.names = F,col.names = F,quote = F)


################################################################################
## Add pseudo color table to result
system(sprintf("(echo %s) | oft-addpct.py %s %s",
               paste0(thres_dir,"color_table.txt"),
               paste0(thres_dir,"tmp_threshold.tif"),
               paste0(thres_dir,"/","tmp_colortable.tif")
))
## Compress final result
system(sprintf("gdal_translate -ot byte -co COMPRESS=LZW %s %s",
               paste0(thres_dir,"/","tmp_colortable.tif"),
               outputfile
))
# gdalinfo(outputfile,hist = T)
## Clean all
system(sprintf(paste0("rm ",thres_dir,"/","tmp*.tif")))

######################################################



#reclassification bfast (Example1)
#bfast_r <- 'example_1.tif'
#outputfile_neg <- 'example_1b2_neg.tif'
outputfile_neg  <- paste0(thres_dir,  substr(basename(bfast_mosaic), 1, nchar(basename(bfast_mosaic))-4),'_nomsk_neg.tif')
#--A_band=2 (no need because the result has already only one band)
system(sprintf("gdal_calc.py -A %s --type=Int16 --NoDataValue=0 --co COMPRESS=LZW --outfile=\"%s\" --calc=%s",
               #paste0(dir,"/",bfast_r),
               result,
               outputfile_neg,
               "\"(A<0)*A\""))


