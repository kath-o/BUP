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
