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

length(unique(longdata$id)) #the number of unique individuals in the dataframe
table(longdata$sex) #equal number of observations of males and females
table(longdata$year) #captures from 1998-2007
table(longdata$island) #at 4 different island locations

#marked requires capture histories to be in 'wide' format
#each row representing one individual
#capture history forms string of 1s and 0s, no spaces, 1 for each capture and 0 where individual isn't detected 
#all histories should have the same number of entrues and therefore be the same length
#should appear in the first column of the dataframe which must be named ‘ch’
#following code uses pipes (%>%) and the spread(), unite() and group_by() functions from the tidyr and dplyr packages to generate the appropriate format

temp <- longdata[,1:2] #take the first two columns, id and year and put into a temporary dataframe
temp$detect <- 1 #add column for detection (all 1s because these represent captures) 

temp <- temp %>%
  #remove duplicates, which may occur when individuals are caught multiple times in an sampling event
  distinct() %>%
  #spread out data. The fill = 0 adds rows for combinations of id and year where individuals were not observed
  spread(year, detect, fill = 0) %>% 
  #for every individual....
  group_by(id) %>%
  #paste together 0's and 1's using unite()
  #here we are pasting the strings together from the second column (first capture event)
  #to the last capture event ("tail(names(.),1)")
  #use sep="" so there are no characters separating 0's and 1's
  unite("ch", 2:tail(names(.),1), sep = "")

sparrow <- as.data.frame(temp) # new dataframe called sparrow
head(sparrow)

#add back some information based on the individual ids using the match() function
sparrow$island <- longdata$island[match(sparrow$id, longdata$id)] 
#this creates a new column called island in the sparrow df...
#using the entry from the island column in the longdata df... 
#where id in the sparrow df matches the id in the longdata df

sparrow$sex <- as.factor(longdata$sex[match(sparrow$id, longdata$id)])

sparrow <- droplevels(subset(sparrow, select = -id)) # remove id column so capture histories appear in first column
head(sparrow)

#Simple Cormack-Jolly-Seber model
#estimates apparent survival (Phi) and detection probability (p) for open populations
#each of these is a linear model on the logit scale (bound between 0 and 1)
#CJS model uses information from capture histories where individuals were marked, not observed, and subsequently seen at a later date to estimate detection probabilities

#most basic form of the model (the default) estimates constant survival and detection probabilities

mod1 <- crm(sparrow) #capture-mark-recapture (cmr) model
mod1 #examine model and coefficient estimates

mod1 <- cjs.hessian(mod1) #refit model with precision estimates

#as with a binomial GLM, these estimates are on the latent (logit) scale. We can transform them back to the data scale using the plogis() or predict() functions
#the estimates on the data scale are also stored within the results section of the model under ‘reals

mod1$results$reals

plogis(mod1$results$beta$Phi)

plogis(mod1$results$beta$p)

predict(mod1, newdata=data.frame(sex = c('Female', 'Male')), se=T) # N.b. In this case, there are no groups or covariates in the model and so the 'newdata' argument is not used

#probability of surviving between capture events was therefore 0.518, and the probability of detecting an animal during a capture event was 0.580

#unequal sampling intervals 
#model assumes an equal time between each capture event
#this assumption can be relaxed by including a vector of time intervals

mod2 <- crm(sparrow, time.intervals = c(1,2,1,1,1,1,1,3,4))
mod2$results$reals

#QUESTIONS: These models assumed constant survival rates and detection probabilities. 
#Is this a realistic assumption for this system? 
#no; both can fluctuate 
#What term might you include in the model next, and why?
#...

#including static covariates
#islands in the meta-population are all different, and it is possible that the probability of capturing an individual differs between islands
#test this by adding some additional complexity to our model and allowing the detection probability to vary between islands

sparrow.proc <- process.data(sparrow) # built in function for data processing
str(sparrow.proc)

head(sparrow.proc[[1]])
head(sparrow.proc$data)

sparrow.ddl <- make.design.data(sparrow.proc) # built in function for building design matrix 
str(sparrow.ddl)
head(sparrow.ddl[[1]])
head(sparrow.ddl$Phi)

#can then specify the model formulation for the detection probability model, making p dependent on island
#we run the CMR model using our newly processed data, our design matrix, and our model specification

# specify model formulation: capture probability depends on island
p.island <- list(formula=~island) 

mod3 <- crm(sparrow.proc, 
            sparrow.ddl, 
            model.parameters = list(p = p.island), 
            accumulate=FALSE, hessian = TRUE)
mod3$results$reals

#Does it look like detection probability varies among islands? 
#Which island has the lowest detection probability? 
#Why might this be? 
#How might you compare this model with our simpler model above to see if it is a better fit?

(mod3$results$AIC)
(mod1$results$AIC)

#lower AIC values indicate a better fitting model, so we accept the model with varying detection probabilities

#test whether survival probabilities differ between the islands
parrow.proc <- process.data(sparrow) 
sparrow.ddl <- make.design.data(sparrow.proc) 

Phi.island <- list(formula=~island) # survival probability depends on island
p.island <- list(formula=~island) # capture probability depends on island

mod4 <- crm(sparrow.proc, 
            sparrow.ddl, 
            model.parameters = list(Phi = Phi.island, 
                                    p = p.island), 
            accumulate=FALSE, hessian = TRUE)

mod4$results$reals
(mod4$results$AIC)
(mod3$results$AIC)

#AIC value of the new model is slightly higher, indicating that including this extra term does not improve model fit. We therefore find no evidence that survival rates vary among islands.

#Does survival probability vary between the sexes?
#test a few different model structures for both survival and detection components
#use a wrapper function 

sparrow.proc <- process.data(sparrow)
sparrow.ddl <- make.design.data(sparrow.proc)

fit.models <- function() {
  Phi.dot <- list(formula=~1) # constant survival
  Phi.sex <- list(formula=~sex) # survival differs between sexes
  Phi.island <- list(formula=~island) # survival differs between islands
  Phi.sex.island <- list(formula=~sex+island) # survival differs between sexes and islands
  p.dot <- list(formula=~1) # constant detection
  p.sex <- list(formula=~sex) # detection probability differs between sexes
  p.island <- list(formula=~island) # detection probability differs between islands
  p.sex.island <- list(formula=~sex+island) # detection probability differs between sexes and islands
  cml <- create.model.list(c("Phi","p"))
  results <- crm.wrapper(cml, data=sparrow.proc, ddl=sparrow.ddl,
                         external=FALSE, accumulate=FALSE, hessian=TRUE)
  return(results)
}

sparrow.models <- fit.models() # run function
sparrow.models # display model table

#top model (with the lowest AIC value) includes a constant model for survival probability with detection probabilities differing among islands
#there are a number of models with very similar AIC scores, indicating a similar level of support
#AIC differences of less than 2 are not considered to be meaningfully different and so in this case we accept the simplest model - that with the fewest parameters

#In this case, that is the top model. We can extract and plot the different detection probabilities as follows.

mod5 <- sparrow.models[[2]]

ggplot(mod5$results$reals$p, aes(island, estimate, ymin=lcl, ymax=ucl)) + 
  geom_errorbar(width=0.2) + geom_point() + ylim(0,1)

#extract the sex and island differences in survival from the similarly supported models and plot those to convince ourselves that there aren’t big differences

mod6 <- sparrow.models[[10]]
ggplot(mod6$results$reals$Phi, aes(sex, estimate, ymin=lcl, ymax=ucl)) + 
  geom_errorbar(width=0.2) + geom_point() + ylim(0,1)

mod7 <- sparrow.models[[6]]
ggplot(mod7$results$reals$Phi, aes(island, estimate, ymin=lcl, ymax=ucl)) + 
  geom_errorbar(width=0.2) + geom_point() + ylim(0,1)

#including time-varying covariates 
#might want to test whether survival (or detection) probabilities differ depending on a factor that varies over time, such as some aspect of the weather
#add a variable ‘cold’ to the design data (sparrow.ddl) and test whether it being an extra-cold winter affects survival

sparrow.ddl$Phi$cold <- "Cold" # new column 
sparrow.ddl$Phi$cold[sparrow.ddl$Phi$time==2 | sparrow.ddl$Phi$time==5 | sparrow.ddl$Phi$time==8] <- "VeryCold" # very cold winters between capture events 2 and 3, 5 and 6, and 8 and 9

head(sparrow.ddl$Phi)

Phi.cold <- list(formula=~cold) 
p.island <- list(formula=~island) 

mod8 <- crm(sparrow.proc, 
            sparrow.ddl, 
            model.parameters = list(Phi = Phi.cold, 
                                    p = p.island), 
            accumulate=FALSE, hessian = TRUE)

mod8$results$reals
(mod8$results$AIC)
(mod5$results$AIC)

ggplot(mod8$results$reals$Phi, aes(cold, estimate, ymin=lcl, ymax=ucl)) + 
  geom_errorbar(width=0.2) + geom_point() + ylim(0,1)

#model including an effect of cold on survival probability has a higher AIC, so there is no evidence that survival varied depending on our made-up variable

#time varying survival and recapture probabilities
# we can test whether survival and detection probabilities varied each year of the study period using the built-in time variable

sparrow.proc <- process.data(sparrow)
sparrow.ddl <- make.design.data(sparrow.proc)

fit.models <- function() {
  Phi.dot <- list(formula=~1) # constant survival
  Phi.time <- list(formula=~time) # survival varies over time
  p.island <- list(formula=~island) # detection probability differs between islands
  p.time <- list(formula=~time) # detection probability varies over time
  p.island.time <- list(formula=~island+time) # detection probability varies over time
  cml <- create.model.list(c("Phi","p"))
  results <- crm.wrapper(cml, data=sparrow.proc, ddl=sparrow.ddl,
                         external=FALSE, accumulate=FALSE, hessian=TRUE)
  return(results)
}

sparrow.models <- fit.models() # run function 

sparrow.models # display model table

mod9 <- sparrow.models[[2]]

(mod9$results$AIC)
(mod5$results$AIC)

ggplot(mod9$results$reals$p, aes(time, estimate, ymin=lcl, ymax=ucl, col=island)) + 
  geom_errorbar(width=0) + geom_point() + ylim(0,1)

#seems that detection probabilities varied over time

#Goodness-of-fit tests
#test whether our data violate any of the basic assumptions of the CJS model, we can use the package R2ucare
#package contains functions to perform goodness-of-fit tests for CMR models
#this package will test whether our data violate the ‘equal detection assumption’ or the ‘equal survival assumption’
#we need to reformat the data into a matrix

sparrow.gof <- sparrow$ch %>%
  strsplit('') %>%
  sapply(`[`) %>%
  t() %>%
  unlist() %>%
  as.numeric %>%
  matrix(nrow = nrow(sparrow))

#R2ucare perfroms three tests, which are unhelpfully called Tests 1, 2 and 3.
#* Test 1: the overall test. Overall, is there evidence that animals have equal detection probabilities and equal survival?
#* Test 2: Does recapture depend on when an animal was first marked? (Tests the equal detection assumption)
#* Test 3: Does marking affect survival? (Tests the equal survival assumption)

overall_CJS(sparrow.gof, rep(1,nrow(sparrow)))

#The p-value is not significant, meaning we fail to reject the null hypothesis. There is therefore no strong evidence for overall lack-of-fit

#If we found evidence of lack-of-fit, we would want to delve deeper to find out what the problem was. First, we can look at the equal detection assumption using Test 2, which has two components:
#* Test 2 CT: Is there a difference in p at t+1 between those captured and not captured at t (when animals are known to be alive because are captured later in the study)?
#* Test 2 CL: Is there a difference in the expected time of next recapture between individuals captured and not captured at t when animals are known to be alive?

test2ct <- test2ct(sparrow.gof, rep(1,nrow(sparrow))) 
test2ct

#Again, we fail to reject the null hypothesis, so no evidence of a problem with the equal detection assumption.

#Finally, Test 3 which tests the equal survival assumption and also has two components:
#* Test 3 SR: Do individuals with previous marks have different survival rates than first-time captures?
#* Test 3 SM: For animals seen again, does when they are recaptured depend on whether they were marked on or before t?

test3sr <- test3sr(sparrow.gof, rep(1,nrow(sparrow)))
test3sr  
  
test3sm <- test3sm(sparrow.gof, rep(1,nrow(sparrow)))
test3sm

#If your data violate the assumptions, you might need to increase the complexity of your model by adding age classes, including a time-varying individual covariate, or building a multistate or multievent model.


