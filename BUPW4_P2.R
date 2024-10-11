#matrix population modelling
#use the popbio package for matrix population modelling

#installing and loading packages 

rm(list=ls(all=TRUE))

my_packages <- c('ggplot2', 'popbio')
new_packages <- my_packages[!(my_packages %in% installed.packages()[,'Package'])]
if(length(new_packages)) install.packages(new_packages, dependencies = T)

library(ggplot2)
library(popbio)


#parameterising your Matrix Population Model 
#stage specific survival 
#CJS model has been run to estimate survival rates of Js, Ys, and As
#comparing a model with and without time-varying detection probabilities, evidence was found that detection probability varied between years
#time-varying detection probabilities to estimate stage-specific survival rates

