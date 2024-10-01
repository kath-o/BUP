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



summary(ouzel_df)