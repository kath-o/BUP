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

mod1 <- crm(sparrow) # capture-mark-recapture (cmr) model
mod1 # examine model and coefficient estimates

mod1 <- cjs.hessian(mod1) # refit model with precision estimates

#As with a binomial GLM, these estimates are on the latent (logit) scale. We can transform them back to the data scale using the plogis() or predict() functions. The estimates on the data scale are also stored within the results section of the model under ‘reals

mod1$results$reals

plogis(mod1$results$beta$Phi)

plogis(mod1$results$beta$p)

predict(mod1, newdata=data.frame(sex = c('Female', 'Male')), se=T) # N.b. In this case, there are no groups or covariates in the model and so the 'newdata' argument is not used

#robability of surviving between capture events was therefore 0.518, and the probability of detecting an animal during a capture event was 0.580

