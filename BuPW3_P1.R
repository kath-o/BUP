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
avi_dat <- read.table('~/Desktop/MSc EEB/WD/BUP/Data/Data_SwissBreedingBirds.csv', header=T, sep=',')

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

#we need to change the projection of our climate data to match that of the bg file.
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

plot(bio_fut)

#3. Model fitting 
#can you code a binomial GLM with all three predictors fitted as linear and squared terms?

model <- glm( Turdus_torquatus ~ bio_2 + I(bio_2^2) + bio_5 + I(bio_5^2) + bio_14 + I(bio_14^2), family='binomial', data=ouzel_df)

summary(model)

#from the model summary, which variables look to be having an important effect?
#bio5 

#4 testing and critiquing our model
#4.1 partial effects 
#model critique is to understand how each climate variable affects the presence/absence of ring ouzels
#we will plot the partial effects, i.e. the effects of each variable whilst holding the other variables constant
#redict function for each variable. In order to use the predict function we will generate a new data frame in which our focal variable varies but the other variables are held constant.

bio_2<-seq(min(ouzel_df$bio_2),max(ouzel_df$bio_2),0.1)

newdata_bio2<-expand.grid(bio_2,mean(ouzel_df$bio_5),mean(ouzel_df$bio_14))
names(newdata_bio2)<-c("bio_2","bio_5","bio_14")
#To use predict we need to generate a new dataset. We will set the other two climate variables at their mean.

response<-predict(model,newdata=newdata_bio2,type="response")
#We've told the predict function to make a prediction on the response scale, ie. in terms of presence/absence

plot(newdata_bio2$bio_2,response,type="l")
