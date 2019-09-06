library(raster)
library(rgdal)
library(sp)

################### List of directories #############
setwd("~")
rootdir <- paste0(getwd(),"/")

aoi_dir <- paste0(rootdir,"Sudan_BFAST/data/AOI_shp/")
bfast_dir <- paste0(rootdir,"Sudan_BFAST/data/bfast_result/")
fo_dir <- ("~/Sudan_BFAST/data/forest_mask_rcl/")



#########################################
#Locate data

aoi<- readOGR(aoi_dir, layer="FREL_STATES")
forestmask <- paste0(fo_dir,'recl_all1234_1302f.tif')   #check raster(...)
bfast_mosaic <- paste0(bfast_dir,'raw_final_results_alltiles.tif')   #check raster(...)
bfast_r <- raster(bfast_mosaic)

crs(raster(forestmask))
crs(raster(bfast_mosaic))

##############################################################################
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

#######################################################################
#reclassify forest mask {
forestmask_rcl <- paste0(fo_dir,'recl_all1234_1302f_wgs_clp.tif')   #check raster(...)

r <- raster(forestmask_rcl)
m <- c(1,11,2, 4,22)
rclmat <- matrix(m, ncol=3, byrow=TRUE)
rc <- reclassify(r, rclmat)

# equivalent to
rc <- reclassify(r, c(1,11,2, 4,22))
# }

plot(rc)

rf <- writeRaster(rc, filename="sdn_forest_shrub_mask_rcl2.tif", format="GTiff", overwrite=TRUE)

forestmask_rcl_sdn <- paste0(fo_dir,'sdn_forest_shrub_mask_rcl2.tif')   #check raster(...)
fmask_rcl <- raster(forestmask_rcl_sdn)

##############################################################################################
# apply forest mask on bfast magnitude product (band 2)

system(sprintf("gdal_calc.py -A %s -B %s --B_band=2 --co COMPRESS=LZW --overwrite --outfile=%s --calc=\"%s\"",
               forestmask_rcl_sdn,
               bfast_mosaic,
               result,
               paste0("A*B")
               
))

result <- raster(result)
plot(result)

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
#divide_sd <- 4
bfast_mg_fmask <- paste0(bfast_dir,"bfast_mg_with_fmask.tif")
result <- raster(bfast_mg_fmask)

tryCatch({
  outputfile <- paste0(bfast_dir,substr(basename(bfast_mosaic), 1, nchar(basename(bfast_mosaic))-4),'_withmsk_threshold.tif')
  means_b2 <- cellStats(result,"mean")
  mins_b2 <- cellStats(result,"min")
  maxs_b2 <- cellStats(result,"max")
  stdevs_b2 <- cellStats(result,"sd")
  system(sprintf("gdal_calc.py -A %s --co=COMPRESS=LZW --type=Byte --outfile=%s --calc='%s'
                   ",
                 bfast_mg_fmask,
                 paste0(bfast_dir,substr(basename(bfast_mosaic), 1, nchar(basename(bfast_mosaic))-4),'_withmsk_threshold.tif'),
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

write.table(pct,paste0(bfast_dir,"color_table.txt"),row.names = F,col.names = F,quote = F)


################################################################################
## Add pseudo color table to result
system(sprintf("(echo %s) | oft-addpct.py %s %s",
               paste0(bfast_dir,"color_table.txt"),
               paste0(bfast_dir,  substr(basename(bfast_mosaic), 1, nchar(basename(bfast_mosaic))-4),'_withfmsk_threshold.tif'),
               paste0(bfast_dir,"/","tmp_colortable.tif")
))
## Compress final result
system(sprintf("gdal_translate -ot byte -co COMPRESS=LZW %s %s",
               paste0(bfast_dir,"/","tmp_colortable.tif"),
               outputfile
))
# gdalinfo(outputfile,hist = T)
## Clean all
system(sprintf(paste0("rm ",bfast_dir,"/","tmp*.tif")))

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




#############################################################
### CLEAN
#system(sprintf("rm %s",
#               paste0(bfast_dir_raw,"tmp_merge.tif")
#))

###################################################################################################################
########### Run change detection
###################################################################################################################


## Perform change detection
system(sprintf("otbcli_MultivariateAlterationDetector -in1 %s -in2 %s -out %s",
               t1_input,
               t2_input,
               imad
))

################################################################################
## Create a no change mask
system(sprintf("gdal_calc.py -A %s  --A_band=4 -B %s  --B_band=1 --outfile=%s --calc=\"%s\"",
               imad,
               imad,
               paste0(imaddir,"tmp_noch.tif"),
               paste0("(A > -0.5)*(A < 0.5)*(B > -0.5)*(B < 0.5)")
))



