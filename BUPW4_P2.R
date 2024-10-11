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
#parameterise a deterministic model
#need a 3x3 matrix with the fertility transitions along the top row, and the survival transitions on the subsequent rows 

# save our estimates of the vital rates
R <- (ClutchNo * HatchingSuc * FledglingNo) / 2
Phi.juv <- survival$estimate[survival$stage=='Juvenile'] 
Phi.yr <- survival$estimate[survival$stage=='Yearling'] 
Phi.ad <- survival$estimate[survival$stage=='Adult'] 

# remind ourselves how these relate to the transition probabilities of the matrix (see slides)
# Juvenile to Juvenile: Phi.juv * R
# Yearling to Juvenile: Phi.yr * R
# Adult to Juvenile: Phi.ad * R
# Juvenile to Yearling: Phi.juv
# Yearling to Yearling: 0 
# Adult to Yearling: 0 
# Juvenile to Adult: 0
# Yearling to Adult: Phi.yr
# Adult to Adult: Phi.ad

# put the transition probabilities into a vector 
sparrowMPM <- c(Phi.juv * R, Phi.yr * R, Phi.ad * R, Phi.juv, 0, 0, 0, Phi.yr, Phi.ad)

# save that vector as a matrix, specifying the number of rows and columns
# use the byrow=TRUE argument to tell R that the first the elements of the vector correspond to the first row of the matrix 
sparrowMPM <- matrix(sparrowMPM, nrow=3, ncol=3, byrow=T)
sparrowMPM

#can now use the popbio package to do some analyses of our deterministic MPM
#for example, we can look at the population growth rate, lambda 

lambda(sparrowMPM) #1.033401

#we can project the dynamics over a given number of iterations (t) based on our matrix and a starting population (n0).

# project over 15 years
t <- 15
# start with 50 juveniles, 20 yearlings and 30 adults
n0 <- c(50,20,30)

# project dynamics 
projection <- pop.projection(sparrowMPM, n0, iterations = t)
projected <- data.frame(time=1:15, N=projection$pop.sizes)

# plot projected pop size over time
ggplot(projected, aes(time, N)) + 
  geom_line() + ylim(0,150) + ylab('Projected N')

#this tells us the population is projected to increase over time since lambda > 1

#observed dynamics 
#we can compare this with estimated population counts (N) over the 15 study years, saved in 'popest.txt'

popest <- read.table("~/Desktop/MSc EEB/WD/BUP/MPM practical/popest.txt", header = TRUE, sep = '\t')
head(popest)

# plot N over time
ggplot(popest, aes(year, N)) + 
  geom_line() + ylim(0,200) + ylab('Observed N')

#stable stage distribution and reproductive value
#popbio has built in functions for other analyses of the asymptotic dynamics of our matrix
#we can look at the stable stage distribution of our population, that is the long term average relative abundance of different stage classes 
#reproductive values of the different stage classes, that is the expected contribution of each individual in that stage class to future reproduction 

stages <- c('Juv','Yr','Ad')
colnames(sparrowMPM) <- stages
rownames(sparrowMPM) <- stages

stable.stage(sparrowMPM)

reproductive.value(sparrowMPM)

#could potentially compare the stable stage distribution with what we observed based on the recapture data to tell us something about the performance of our model 

#perturbation analysis 
#popbio also includes built-in functions to calculate the sensitivities and elasticities of the different vital rates 
#these tell us about the relative importance of each vital rate in determining the population growth rate, lambda
#sensitivities estimate the change in lambda for an absolute change in a vital rate 
#elasticities tell us about the effect of a proportional change 

# list the vital rates
sparrow.param <- list(Phi.juv = Phi.juv, Phi.yr = Phi.yr, Phi.ad = Phi.ad, R = R)

# give the matrix equation 
sparrow.equation <- expression(Phi.juv * R, Phi.yr * R, Phi.ad * R, Phi.juv, 0, 0, 0, Phi.yr, Phi.ad)

# run the sensitivity analysis
sens <- vitalsens(sparrow.equation, sparrow.param)
sens

# plot elasticity of the vital rates 
sens$vitalrate <- factor(c('Phi.juv', 'Phi.yr', 'Phi.ad', 'R'), levels = c('Phi.juv', 'Phi.yr', 'Phi.ad', 'R'))
ggplot(sens, aes(vitalrate, elasticity)) + 
  geom_bar(stat = 'identity') 


