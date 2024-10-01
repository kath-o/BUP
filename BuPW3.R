#new project 
#sorting git 

#installing raster
install.packages("raster",dependencies=TRUE,repos="https://cloud.r-project.org")
library(raster)

#loading data
avi_dat <- read.table('Data_SwissBreedingBirds.csv', header=T, sep=',')

nrow(avi_dat)

summary(avi_dat)

#subset the data to select columns were working with 

ouzel_cols <- c('Turdus_torquatus', 'bio_5', 'bio_2', 'bio_14', 'blockCV_tile')

ouzel_df <- data.frame(avi_dat)[ouzel_cols]

summary(ouzel_df)
