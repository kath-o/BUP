#installing packages 
install.packages("geodata",dependencies=TRUE,repos="https://cloud.r-project.org")

install.packages("predicts",dependencies=TRUE,repos="https://cloud.r-project.org")

install.packages("terra",dependencies=TRUE,repos="https://cloud.r-project.org")

library(geodata)
library(predicts)
library(terra)

#downloading data from GBIF - first 10000, my species rimula latifolia
occdata <- geodata::sp_occurrence("Primula", "latifolia*", geo=FALSE,removeZeros=TRUE,start=1,end=10000)

#star gives us all name variants
#geo = TRUE only shows records w lat and long 

#data frame
dim(occdata)
occdata[1:10, ]

#plot global distribution to make sure it fits with expectation 
wrld <- world(path="~/Desktop/MSc EEB/WD/BUP")
#this function gives us an outline of the world's political boundaries. Reminder, if ever you want to know more about an R function, you can write ?function.name, e.g., ?world
plot(wrld, xlim=c(-180,180), ylim=c(-80,80), col="light yellow", border="light gray")
# add the points
points(occdata$lon, occdata$lat, col='blue', pch=20)

#cleaning up occurrence data 
#GBIF can have errors, needs to be scrutinised 
#not necessary for this data, but if it was: 
#occdata<-subset(occdata,lat>0)
#dups <- duplicated(occdata[, c('lon', 'lat')])
#This identifies observations that have already appeared above
#sum(dups)
#There are a lot of them, so removing them will give us a smaller dataset to work with
#occ <- occdata[!dups, ]

#downloading worldclim data 
#download data at a resolution of 10 minutes (60 mins a degree)

output_dir<-"~/Desktop/MSc EEB/WD/BUP/Data"

bio_glob<-worldclim_global(var="bio", res=10,path=output_dir, version="2.1")

dim(bio_glob)

#bio_glob is in the form of a spatraster (i.e. spatial raster). A raster is simply another name for a spatial grid, and a spatraster combines several of these grids - one for each of the 19 worldclim variables
#we will also clip the spatraster so it only covers the spatial extent of our study species. First its longitudes then latitudes
summary(occ$lon)
summary(occ$lat)

e <- ext(-5, 20, 30, 60)

predictors <- crop(bio_glob, e)

names(predictors)<-substring(names(predictors),11,16)

#we can now have a look at the global climate data. Here we’ll just look at the first 9 worldclim variables
plot(predictors,1:9)

#and here we can add our species data onto a plot of climate data for the first variable.
plot(predictors,19)
points(occ$lon,occ$lat, col='maroon2',pch=16,cex=0.2)

#generating background data
#occurrence only data: sample background data from a region, should cover locations where the species is present and absence

#here I'm setting the spatial extent to be broadly consistent with that of my study species (you need to make sure it is sampling from the same extent). Remember to find out how a function works you can do ?function
bg<-spatSample(predictors,5000,"random", na.rm=TRUE, as.points=TRUE,ext=e)

#Here we'll plot our background points on a map of climwin variable 1 (you could change this to any of the worldclim variables)
plot(predictors, 1)
points(bg, cex=0.1, col="darkmagenta")

#matching occurrence and climate data
occlatlon<-cbind(occ$lon,occ$lat)
presvals <- extract(predictors, occlatlon)
#presvals is the climate data for where the species is present
backvals <- values(bg)
#backvals is the climate data for the background data
bg_lonlat<-geom(bg)
lonlats<-rbind(occlatlon, bg_lonlat[,c("x","y")])
pb <- c(rep(1, nrow(presvals)), rep(0, nrow(backvals)))
#The first column of the dataset is a vector of 1s for presences and 0s for background data.
sdmdata <- data.frame(cbind(lonlats,pb, rbind(presvals, backvals)))
#here we combine the presence and background data into a single data frame

#in sdmdata the first two columns are climate data and third column contains either a 1 (species present) or 0 (background data). The remaining columns are the corresponding worldclim variables for these location.

#we can also examine how colinear (i.e. correlated) predictor variables are. Highly correlated predictor variables can give rise to statistical issues.
pairs(sdmdata[,4:7], cex=0.1)
#here we just look at the correlations between the first 4 climate variables. You could extend this to look at all 19. 

#fitting a species distribution model 
#the predicts package in R allows us to fit a whole variety of functions that are useful in running species distribution models. Next we will use it to run the MAXENT approach that has proven very popular.
sdmdata<-subset(sdmdata,is.na(bio_1)==F)
#here we're just removing a couple of rows where the climate data are NAs.

specdata<-as.data.frame(cbind(rep("Primula latifola",length(sdmdata[,1])),
                              sdmdata))

names(specdata)[1:4]<-c("species","longitude","latitude","presabs")

specdata<-subset(specdata,presabs==1)

backdata<-as.data.frame(cbind(rep("background",length(sdmdata[,1])),
                              sdmdata))

names(backdata)[1:4]<-c("","longitude","latitude","presabs")

backdata<-subset(backdata,presabs==0)


write.table(specdata[,-4],paste(output_dir,"/Primulalatifola.csv",sep=""),col.names=T,row.names=F,sep=",")
write.table(backdata[,-4],paste(output_dir,"/background.csv",sep=""),col.names=T,row.names=F,sep=",")

model<-MaxEnt(sdmdata[,-c(1:3)],sdmdata[,3],removeDuplicates=TRUE)

#here we've used all of the climate variables, but you could be more discerning. 
#we've also asked the model to ignore any data that comes from the same cell.
#in the maxent call we first specify the climate data - these are in all the columns except the first one.
#next we specify the presence/background data - column1 i.e. [,3]

plot(model)
model

#We can look at the predicted climate suitability globally, and see how it matches where the species has been recorded (remember with occurrence data no data doesn’t mean it is necessarily absent).
predictedocc <- predict(model, predictors, args=c("outputformat=raw")) 

par(mfrow=c(2,1))
plot(predictedocc)
plot(predictedocc)
points(occlatlon,pch=".", col="deeppink")

#predicting future distributions 
#we need to download future climate data

bio_fut<-cmip6_world(model='ACCESS-ESM1-5', ssp='245', time='2041-2060', var='bioc', res=10, path=output_dir)

fut_predictors<-crop(bio_fut,e)

plot(predictors,17)

plot(fut_predictors,17)

#generate a prediction of future climate suitability for our species
names(fut_predictors)<-names(predictors)

fut_predictedocc <- predict(model, fut_predictors, args=c("outputformat=raw")) 

par(mfrow=c(2,1))
plot(predictedocc,main="current")

plot(fut_predictedocc,main="2050")
