####################################################################################################
####################################################################################################
## Tiling of an AOI (shapefile defined)
## Contact remi.dannunzio@fao.org 
## 2019/03/11
####################################################################################################
####################################################################################################
############ CREATE A FUNCTION TO GENERATE REGULAR GRIDS
generate_grid <- function(aoi_utm,size){
  ### Create a set of regular SpatialPoints on the extent of the created polygons  
  sqr <- SpatialPoints(makegrid(aoi_utm,offset=c(-0.5,-0.5),cellsize = size))
  
  ### Convert points to a square grid
  grid <- points2grid(sqr)
  
  ### Convert the grid to SpatialPolygonDataFrame
  SpP_grd <- as.SpatialPolygons.GridTopology(grid)
  
  sqr_df <- SpatialPolygonsDataFrame(Sr=SpP_grd,
                                     data=data.frame(rep(1,length(SpP_grd))),
                                     match.ID=F)
  ### Assign the right projection
  proj4string(sqr_df) <- proj4string(aoi)
  sqr_df
}

#plot(sqr_df, add=TRUE)
### load the parameters
#source('~/uga_activity_data/scripts/get_parameters.R')
#usernamelist <- paste0(mgmt_dir,'usernames_uga.csv')


### GET COUNTRY BOUNDARIES FROM THE WWW.GADM.ORG DATASET
#aoi   <- getData('GADM',
#                 path=gadm_dir, 
#                 country= countrycode, 
#                 level=0)

#(bb    <- extent(aoi))
aoi <- readOGR(dsn="~/4_sudan/data/aoi/Sudan_States.shp")
crs(aoi)

plot(aoi)
aoiutm <- spTransform(aoi, CRS("+proj=utm +zone=36 ellps=WGS84"))
aoi <- spTransform(aoi, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))
aoiutm

 
### What grid size do we need ? 
#grid_size <- 20000          ## in meters
grid_size <- 120000           ## in meters (36 tiles)
### GENERATE A GRID 
sqr_df <- generate_grid(aoi,grid_size/111320) ## in degree
aoi@data
nrow(sqr_df)
plot(sqr_df)
### Select a vector from location of another vector
sqr_df_selected <- sqr_df[aoi,]
nrow(sqr_df_selected)

### Give the output a decent name, with unique ID
#names(sqr_df_selected@data) <- "tileID" 
names(sqr_df_selected@data) <- "tileIDm" 
head(sqr_df_selected)
#sqr_df_selected@data$tileID <- row(sqr_df_selected@data)[,1]
sqr_df_selected@data$tileIDm <- row(sqr_df_selected@data)[,1]

tiles   <- spTransform(sqr_df_selected,CRS("+init=epsg:4326"))
#aoi_geo <- spTransform(aoi,CRS("+init=epsg:4326"))

#include states (ensure stame crs)
tiles <- spTransform(tiles,CRS("+proj=longlat +datum=WGS84"))
aoi <- spTransform(aoi,CRS("+proj=longlat +datum=WGS84"))
tiles@data$state <- over(tiles,aoi)$state

#tilesstates <- SpatialPointsDataFrame(c...
head(tiles)
#rename cols

#names(tiles@data) <- c('tileID','State')
names(tiles@data) <- c('tileID','tileIDm')
head(tiles@data)


tiles(data)[tiles(data)=="State.State"] <- "State"
# Rename a column in R
names(tiles)[2]<-"State"

### Reproject in LAT LON
tiles   <- spTransform(sqr_df_selected,CRS("+init=epsg:4326"))
aoi_geo <- spTransform(aoi,CRS("+init=epsg:4326"))


### Plot the results
plot(tiles, add=T)
plot(aoi_geo,add=T,border="blue")


### check against forest nonforest mask
#lc <- paste0(lc_dir,'recl_all1234_1302f.tif')
#lc_wgs <- spTransform(lc, CRS("+proj=utm +zone=36 ellps=WGS84"))
#aoi <- spTransform(aoi, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))

FNF_mask_proj <- paste0(lc_dir,"recl_all1234_1302frcl.tif")
FNF_mask_proj.shp <- paste0(lc_dir,"recl_all1234_1302frcl.shp")
shp_mask <- readOGR(dsn="~/4_sudan/data/forest_mask/recl_all1234_1302frcl.shp")
#aoi <- readOGR(dsn="~/4_sudan/data/aoi/Sudan_States.shp")
#plot(FNF_mask_proj.shp, add=T, border="red")

#FNF_mask_proj <- paste0(lc_dir,"FNF_mask_2015_2017_proj.tif")
#FNF_mask_proj.shp <- paste0(lc_dir,"FNF_mask_2015_2017_proj.shp")
plot(raster(FNF_mask_proj),add=T)
######################################
if(!file.exists(FNF_mask_proj.shp)){
system(sprintf("gdal_polygonize.py -mask %s %s -f 'ESRI Shapefile' %s  ",
                               FNF_mask_proj,
                               FNF_mask_proj,
                               FNF_mask_proj.shp
                               ))
}

#shp_mask <- readOGR(FNF_mask_proj.shp)
# plot(shp_mask,add=T)
shp_mask <- spTransform(shp_mask,CRS("+init=epsg:4326"))

tiles@data$forest_mask <- over(tiles,shp_mask)
subtile <- tiles[tiles@data$forest_mask$DN %in% 1 ,]
subtile <- subtile[,"tileID"]
subtile
table(tiles$forest_mask)

subtile <- tiles[tiles@data$forest_mask$DN %in% 1 ,]
subtile <- subtile[,"tileID"]
plot(shp_mask,border="green")
plot(subtile,add=T, border="red")
plot(subtile)

### Read the list of usernames
users     <- read.csv(usernamelist)

### Assign each tile with a username
df        <- data.frame(cbind(subtile@data[,"tileID"],users$Name))
names(df) <- c("tileID","username")
df$tileID <- as.numeric(df$tileID)
table(df$username)
subtile@data <- df


for(username in unique(df$username)){
  ### Create a final subset corresponding to your username
    my_tiles <- subtile[subtile$tileID %in% df[df$username == username,"tileID"],]
    # plot(my_tiles,add=T,col="black")
    length(my_tiles)
    
    ### Export the final subset
    export_name <- str_replace_all(paste0("national_scale_",length(my_tiles),"_tiles_",username), " ","_")
    
    
    writeOGR(obj=my_tiles,
        dsn=paste(tile_dir,export_name,".kml",sep=""),
        layer= export_name,
        driver = "KML",
        overwrite_layer = T)
    }
    
### Export ALL TILES as KML
#export_name <- paste0("tiling_system_",countrycode)
#export_name <- "tiling_system_SDN"
export_name <- "entire_country"

#writeOGR(obj=my_tiles,
#         dsn=paste(tile_dir,export_name,".kml",sep=""),
#         layer= export_name,
#         driver = "KML",
#         overwrite_layer = T)

writeOGR(obj=subtile,
         dsn=paste(tile_dir,export_name,".kml",sep=""),
         layer= export_name,
         driver = "KML",
         overwrite_layer = T)
#re-run gws transfor of the tiles var
writeOGR(obj=tiles,
         dsn=paste(tile_dir,export_name,".kml",sep=""),
         layer= export_name,
         driver = "KML",
         overwrite_layer = T)


my_tiles <- tiles[tiles$tileID %in% df[df$username == username,"tileID"],]

plot(my_tiles,add=T,col="green")
length(my_tiles)


### Export the ONE TILE IN THE subset
export_name <- paste0("UGA_one_tile_",username)

writeOGR(obj=my_tiles[1,],
         dsn=paste(tile_dir,export_name,".kml",sep=""),
         layer= export_name,
         driver = "KML",
         overwrite_layer = T)

#subset 10 tiles randomly
n10_tiles = subtile[sample(2431,10),]
plot(aoi)
plot(n10_tiles, add=TRUE)
export_name <- paste0("tiling_system_6v2x_",countrycode)
writeOGR(obj=n10_tiles,
         dsn=paste(tile_dir,export_name,".kml",sep=""),
         layer= export_name,
         driver = "KML",
         overwrite_layer = T)

writeOGR(obj=tiles,
         dsn=paste(tile_dir,export_name,".kml",sep=""),
         layer= export_name,
         driver = "KML",
         overwrite_layer = T)


###################
##merge the two shapefiles (tile 20 km and tile 60 km )

library(rgeos)
library(maptools)
### Plot the results
plot(tiles, add=T)
plot(aoi_geo,add=T,border="blue")

#tile_clp_country = gIntersection(aoi_geo,tiles)
#plot(tile_clp_country )

rgdal::ogrListLayers("foo.kml")
rootdir <- paste0(getwd(),"/")
bfast_tile_dir<- paste0(rootdir,"data/bfast/tiles/")


#tile20kminfo <- system(paste("ogrinfo", paste(bfast_tile_dir, "tiling_system_SDN.kml", sep="/")), intern=TRUE)
#tiles20km <- paste0(bfast_tile_dir,'tiling_system_SDN.kml')
tiles20km  <- readOGR(paste(bfast_tile_dir, "tiling_system_SDN.kml", sep="/"), "tiling_system_SDN")
bb <- bbox(tiles20km)
bb
plot(tiles20km)
e <- extent(aoiutm)
# Create extent and coerce to SpatialPolygons 
# Create extent and coerce to SpatialPolygons 
ee <- as( raster::extent(-721125.8, 1091414, 966396.9, 2562199), "SpatialPolygons")

#e <- as( raster::extent(334367, 334498.7, 4088915, 4089057), "SpatialPolygons")
proj4string(ee) <- "+proj=utm +zone=36 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
#coord. ref. : +proj=utm +zone=36 ellps=WGS84 +ellps=WGS84 
class(e)
plot(ee)

##############################################################
# Create raster from defined extent
#https://gis.stackexchange.com/questions/261054/setting-a-grid-inside-a-shapefile-in-r
library(raster)
library(sp)
e
aoi
cs= 20000
r <- raster(ee, resolution = cs)
r[] <- 1:ncell(r)
plot(r)
cat("\n", "Number of cells in raster: ", ncell(r), "\n")  
# Coerce to SpatialPixelsDataFrame
r <- as(r, "SpatialPixelsDataFrame")
class(r)
# Plot cell boundaries (grid) and SpatialPixelsDataFrame raster
plot( as(r, "SpatialPolygons") )
plot(r)
plot(aoiutm, add=TRUE)

cs2= 120000
r_cs2 <- raster(ee, resolution = cs2)
r_cs2[] <- 1:ncell(r_cs2)
plot(r_cs2, add=TRUE)
cat("\n", "Number of cells in raster: ", ncell(r_cs2), "\n")  
r_cs2 <- as(r_cs2, "SpatialPixelsDataFrame")
class(r_cs2)
plot( as(r_cs2, "SpatialPolygons") )
plot( as(r, "SpatialPolygons"), add=TRUE)



### Give the output a decent name, with unique ID
#names(sqr_df_selected@data) <- "tileID" 
names(r@data) <- "tileID20" 
names(r_cs2@data) <- "tileID120"
#sqr_df_selected@data$tileID <- row(sqr_df_selected@data)[,1]
r@data$tileID20 <- row(r@data)[,1]
r_cs2@data$tileID120 <- row(r_cs2@data)[,1]
#https://gis.stackexchange.com/questions/155328/merging-multiple-spatialpolygondataframes-into-1-spdf-in-r
r_r_cs2 <-union(r, r_cs2)
joined = SpatialPolygons(lapply(spl, function(x){x@polygons[[1]]}))
head(r_r_cs2)
ab<-merge(a, b)
ab <- bind(r,r_cs2,makeUniqueIDs = TRUE)
head(ab)
plot(ab)
export_name <- paste0("ab_6v2x_",countrycode)


abwgs <- spTransform(ab,CRS("+proj=longlat +datum=WGS84"))

writeOGR(obj=abwgs,
         dsn=paste(tile_dir,export_name,".kml",sep=""),
         layer= export_name,
         driver = "KML",
         overwrite_layer = T)

b <- spChFIDs(r, paste("99", row.names(r), sep="."))


#https://gis.stackexchange.com/questions/180682/merge-a-list-of-spatial-polygon-objects-in-r?rq=1
#http://www.petrkeil.com/?p=648
#http://www.nickeubank.com/wp-content/uploads/2015/10/RGIS2_MergingSpatialData_part1_Joins.html


r@data[c("tileID20")]

bo <- merge(r, r_cs2, by.x = "tileID20", by.y = "tileID120", duplicateGeoms = TRUE)   #if TRUE geometries in x are duplicated if there are multiple matches between records in x and y
duplicated(r_cs2$tileID120)
duplicated(r$tileID20)
bo <- merge(r_cs2,r, by.x = "tileID120", by.y = "tileID20", duplicateGeoms = TRUE)   #if TRUE geometries in x are duplicated if there are multiple matches between records in x and y
#https://stackoverflow.com/questions/38256946/what-can-cause-a-non-unique-matches-detected-error-in-an-r-merge

tiles60km  <- readOGR(paste(bfast_tile_dir, "tiling_system_6v2x_SDN.kml", sep="/"), "tiling_system_6v2x_SDN")    #tiling_system_6x_SDN
plot(tiles20km)
subs_union <- union(subs1, subs2)
xxx <- union(tiles20km, tiles60km)
xx <- bind(tiles20km, tiles60km)
head(xx)


plot(x)
head(tiles)
head(tiles20km)
head(tiles60km)
export_name <- paste0("tiling_sys62_final_",countrycode)
writeOGR(obj=xx,
         dsn=paste(tile_dir,export_name,".shp",sep=""),
         layer= export_name,
         driver = "KML",
         overwrite_layer = T)

