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
#geo = TRUE only shows 
