####################################################################################################
####################################################################################################
## Tiling of an AOI (shapefile defined)
## Contact remi.dannunzio@fao.org 
## 2019/03/11
####################################################################################################
####################################################################################################
############ CREATE A FUNCTION TO GENERATE REGULAR GRIDS
library(rgdal)
library(stringr)
library(recipes)

lc_dir <- "~/downloads/Sudan_BFAST/Sudan/"
rootdir       <- "~/downloads/Sudan_BFAST/Sudan/"
tile_dir <- "~/downloads/Sudan_BFAST/Sudan/"
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
mgmt_dir <-('~/downloads/Sudan_BFAST/Sudan/')
usernamelist <- paste0(mgmt_dir,'participants_workshop_general.csv')


### GET COUNTRY BOUNDARIES FROM THE WWW.GADM.ORG DATASET

aoi <- readOGR(dsn="~/downloads/Sudan_BFAST/Sudan/Sudan_States.shp")

crs(aoi)
plot(aoi)

aoiutm <- spTransform(aoi, CRS("+proj=utm +zone=36 ellps=WGS84"))
aoi <- spTransform(aoi, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))

 
### What grid size do we need ? 
#grid_size <- 20000          ## in meters
grid_size <- 40000           ## in meters (40 x 40)
### GENERATE A GRID 
sqr_df <- generate_grid(aoi,grid_size/111320) ## in degree
aoi@data
nrow(sqr_df)
plot(sqr_df)
### Select a vector from location of another vector
sqr_df_selected <- sqr_df[aoi,]

nrow(sqr_df_selected)

### Give the output a decent name, with unique ID

names(sqr_df_selected@data) <- "tileIDmm" 
head(sqr_df_selected)

sqr_df_selected@data$tileIDmm <- row(sqr_df_selected@data)[,1]

tiles   <- spTransform(sqr_df_selected,CRS("+init=epsg:4326"))


#include states (ensure same crs)
tiles <- spTransform(tiles,CRS("+proj=longlat +datum=WGS84"))
aoi <- spTransform(aoi,CRS("+proj=longlat +datum=WGS84"))
tiles@data$state <- over(tiles,aoi)$state


#tilesstates <- SpatialPointsDataFrame(c...
head(tiles)
#rename cols

#names(tiles@data) <- c('tileIDm','State')
#names(tiles@data) <- c('tileIDm','tileIDm')
head(tiles@data)


# tiles(data)[tiles(data)=="State.State"] <- "State"
# # Rename a column in R
# names(tiles)[2]<-"State"

### Reproject in LAT LON
tiles   <- spTransform(sqr_df_selected,CRS("+init=epsg:4326"))
aoi_geo <- spTransform(aoi,CRS("+init=epsg:4326"))


### Plot the results
plot(tiles, add=T)
plot(aoi_geo,add=T,border="blue")


### check against forest nonforest mask

FNF_mask_proj <- paste0(lc_dir,"recl_all1234_1302frcl.tif")
FNF_mask_proj.shp <- paste0(lc_dir,"recl_all1234_1302frcl.shp")
shp_mask <- readOGR(dsn="~/downloads/Sudan_BFAST/Sudan/recl_all1234_1302frcl.shp")
FRL_states <- readOGR(dsn="~/downloads/Sudan_BFAST/Sudan/FREL_STATES.shp")
FRL_states <- spTransform(FRL_states,CRS("+init=epsg:4326"))
#plot(FNF_mask_proj.shp, add=T, border="red")

######################################
if(!file.exists(FNF_mask_proj.shp)){
system(sprintf("gdal_polygonize.py -mask %s %s -f 'ESRI Shapefile' %s  ",
                               FNF_mask_proj,
                               FNF_mask_proj,
                               FNF_mask_proj.shp
                               ))
}

shp_mask <- readOGR(FNF_mask_proj.shp)
plot(shp_mask,add=T)
shp_mask <- spTransform(shp_mask,CRS("+init=epsg:4326"))

tiles@data$forest_mask <- over(tiles,shp_mask)
subtile <- tiles[tiles@data$forest_mask$DN %in% 1 ,]
subtile <- subtile[,"tileIDmm"]
subtile
table(tiles$forest_mask)

#subtile <- tiles[tiles@data$forest_mask$DN %in% 1 ,]
#subtile <- subtile[,"tileIDm"]
plot(shp_mask,border="green")
plot(subtile,add=T, border="red")
plot(subtile)

####Clip to FRL states only

tiles_FRL_states <- subtile[FRL_states, ]
plot(tiles_FRL_states)

### Read the list of usernames
users     <- read.csv(usernamelist)

### Assign each tile with a username
df        <- data.frame(cbind(tiles_FRL_states@data[,"tileID"],users$Name)) # number of rows of result is not a multiple of vector length (arg 2)
names(df) <- c("tileID","username")
df$tileID <- as.numeric(df$tileID)
table(df$username)


nbatch <- ceiling(max(table(df$username)))

df <- arrange(df,username)
df <- cbind(df,rep(1:nbatch,ceiling(nrow(df)/nbatch))[1:nrow(df)])
names(df) <- c("tileID","username","batch")
df <- arrange(df,tileID)

#table(df$username,df$batch)
tiles_FRL_states@data <- df


for(i in 1:dim(users)[1]){
  #for(username in users){
  #for(batch in nbatch){
  ### Create a final subset corresponding to your username
  my_tiles <- tiles_FRL_states[tiles_FRL_states$tileID %in% df[df$username == users[i,],"tileID"], ]
  #plot(my_tiles,add=T,col="red")
  length(my_tiles)
  
  ### Export the final subset
  
  export_name <- str_replace_all(paste0("national_scale_",length(my_tiles),"_tiles_",users[i,]), "batch ","_")
  
  rgdal:: writeOGR(obj=my_tiles,
                   dsn=paste(tile_dir,export_name,".kml",sep=""),
                   layer= export_name,
                   driver = "KML",
                   overwrite_layer = F)
  
}


   
