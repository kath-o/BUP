#new project 
#sorting git 
#fixed error using "git config --global http.postBuffer 524288000" - ask

#installing raster
install.packages("raster",dependencies=TRUE,repos="https://cloud.r-project.org")
library(raster)

#change to geodata and terra

install.packages("geodata",dependencies=TRUE,repos="https://cloud.r-project.org")
library(geodata)
library(terra)

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
