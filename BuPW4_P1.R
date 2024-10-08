#packages 
my_packages <- c('dplyr', 'tidyr', 'marked', 'ggplot2', 'R2ucare')
new_packages <- my_packages[!(my_packages %in% installed.packages()[,'Package'])]
if(length(new_packages)) install.packages(new_packages)

library(dplyr)
library(tidyr)
library(marked)
library(ggplot2)
library(R2ucare)

#load the sparrow recapture dataset
longdata <- read.table("~/Desktop/MSc EEB/WD/BUP/CMR practical/sparrowrecap.txt", header = TRUE, sep = '\t')
head(longdata)

#dataset consists of n=2378 capture observations across different islands in the house sparrow meta-population
#ata are in what is known as ‘long format’, with individuals appearing in multiple rows representing each occasion they were captured
#simple data exploration

length(unique(longdata$id)) # the number of unique individuals in the dataframe
table(longdata$sex) # equal number of observations of males and females
table(longdata$year) # captures from 1998-2007
table(longdata$island) # at 4 different island locations
