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

#Juvenile survival: 0.463 (95% CI 0.404–0.524)
#Yearling survival: 0.510 (95% CI 0.445–0.574)
#Adult 2+ survival: 0.559 (95% CI 0.499–0.618)

#plotting these values 

# first enter the values into a dataframe 
survival <- data.frame(stage=factor(c('Juvenile','Yearling','Adult'), levels=c('Juvenile','Yearling','Adult')), estimate=c(0.463, 0.510, 0.559), lcl=c(0.404, 0.445, 0.499), ucl=c(0.524, 0.574, 0.618))

# then plot by stage
ggplot(survival, aes(stage, estimate, ymin=lcl, ymax=ucl)) + 
  geom_errorbar(width=0.2) + geom_point() + ylim(0,1)

#juveniles seem to have the lowest survival rate, which is to be expected 

#per capita reproduction 
nestdata <- read.table("~/Desktop/MSc EEB/WD/BUP/MPM practical/gjeroynest.txt", header = TRUE, sep = '\t')
head(nestdata)

#clutchno indicates whether it was the first, second, or third clutch laid in that nest in that breeding season
#hatchingsuc 

#estimating per capita reproduction from these data 
HatchingSuc <- mean(nestdata$hatchingsuc)
FledglingNo <- mean(nestdata$chickno)

#for the number of clutches, create a new data frame called nest; one row for each unique nest
#add a column for this, which takes the maximum value of clutchno for each unique value of nestid 
#we then take the mean of these values to be the average number of clutches 

nests <- data.frame(nestid = sort(unique(nestdata$nestid)), numberofclutches=tapply(nestdata$clutchno, nestdata$nestid, max))
ClutchNo <- mean(nests$numberofclutches)

#we can use these as estimates of clutch number, hatching success and fledgling number to claculate per capita reproduction
#then take average number of clutches (clutch number) and multiply by the probability of a clutch hatching (hutching succ), multiplied by the expected number of fledglings (fledgling number)
#gives us the expected number of chicks per female over the breeding season 
#since we are modelling only the female segment of the population, we will then divide this estimate by 2 to get the number of female chicks, assuming an equal sex ratio of offspring 

(ClutchNo * HatchingSuc * FledglingNo) / 2

#deterministic population model 




