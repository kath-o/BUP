#new project 
#sorting git 
#fixed error using "git config --global http.postBuffer 524288000" - ask

#1. Conceptualisation
#installing raster
install.packages("raster",dependencies=TRUE,repos="https://cloud.r-project.org")
library(raster)

#change to geodata and terra

install.packages("geodata",dependencies=TRUE,repos="https://cloud.r-project.org")
library(geodata)
library(terra)

#2. Data preparation
#loading data
avi_dat <- read.table('Data_SwissBreedingBirds.csv', header=T, sep=',')

nrow(avi_dat)

summary(avi_dat)

#subset the data to select columns were working with 

ouzel_cols <- c('Turdus_torquatus', 'bio_5', 'bio_2', 'bio_14', 'blockCV_tile')

ouzel_df <- data.frame(avi_dat)[ouzel_cols]

summary(ouzel_df)

#importing current and future climate data for Switzerland and clip it to the extent of Switzerland
#Download and import same climate variables from WORLDCLIM using geodata package

output_dir<-"~/Desktop/MSc EEB/WD/BUP/Data"

bio_curr <-worldclim_country("Switzerland",version="2.1", var='bio', res=10, lon=5.5, lat=45.5, path=output_dir)[[c(2,5,14)]]

bio_fut <- cmip6_world(var = "bio", model = "CNRM-CM6-1-HR", ssp = "245", res = 10,  time = "2041-2060",  lon = c(5.96, 10.49),  lat = c(45.82, 47.81),path=output_dir)[[c(2,5,14)]]

#next step is to reproject the raster data in Swiss coordinates and first crop the layers for the Swiss boundary

# A spatial mask of Switzerland in Swiss coordinates
bg <- rast('/vsicurl/https://damariszurell.github.io/SDM-Intro/CH_mask.tif')

bio_curr <- terra::project(bio_curr, bg)
bio_fut <- terra::project(bio_fut, bg)
#we need to change the projection of our cliamte data to match that of the bg file.

bio_curr <- terra::resample(bio_curr, bg)
bio_fut <- terra::resample(bio_fut, bg)
#we then need to make the resolution equivalent to bg. 


bio_curr <- terra::mask(bio_curr, bg)
bio_fut <- terra::mask(bio_fut, bg)
#we then need to clip the extent to match an outline of Switzerland

names(bio_curr) <- c('bio_2', 'bio_5', 'bio_14')
names(bio_fut) <- c('bio_2', 'bio_5', 'bio_14')

#plot the current climate variables and how they are projected to change over time
plot(bio_curr)
