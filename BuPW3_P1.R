#new project 
#sorting git 
#fixed error using "git config --global http.postBuffer 524288000" - ask

#1. Conceptualisation
#installing raster
install.packages("raster",dependencies=TRUE,repos="https://cloud.r-project.org")
library(raster)

#change to geodata and terra

install.packages("geodata",dependencies=TRUE,repos="https://cloud.r-project.org")
library(geodata)
library(terra)

#2. Data preparation
#loading data
avi_dat <- read.table('~/Desktop/MSc EEB/WD/BUP/Data/Data_SwissBreedingBirds.csv', header=T, sep=',')

nrow(avi_dat)

summary(avi_dat)

#subset the data to select columns were working with 

ouzel_cols <- c('Turdus_torquatus', 'bio_5', 'bio_2', 'bio_14', 'blockCV_tile')

ouzel_df <- data.frame(avi_dat)[ouzel_cols]

summary(ouzel_df)

#importing current and future climate data for Switzerland and clip it to the extent of Switzerland
#Download and import same climate variables from WORLDCLIM using geodata package

output_dir<-"~/Desktop/MSc EEB/WD/BUP/Data"

bio_curr <-worldclim_country("Switzerland",version="2.1", var='bio', res=10, lon=5.5, lat=45.5, path=output_dir)[[c(2,5,14)]]

bio_fut <- cmip6_world(var = "bio", model = "CNRM-CM6-1-HR", ssp = "245", res = 10,  time = "2041-2060",  lon = c(5.96, 10.49),  lat = c(45.82, 47.81),path=output_dir)[[c(2,5,14)]]

#next step is to reproject the raster data in Swiss coordinates and first crop the layers for the Swiss boundary

# A spatial mask of Switzerland in Swiss coordinates
bg <- rast('/vsicurl/https://damariszurell.github.io/SDM-Intro/CH_mask.tif')

bio_curr <- terra::project(bio_curr, bg)
bio_fut <- terra::project(bio_fut, bg)

#we need to change the projection of our climate data to match that of the bg file.
bio_curr <- terra::resample(bio_curr, bg)
bio_fut <- terra::resample(bio_fut, bg)

#we then need to make the resolution equivalent to bg. 
bio_curr <- terra::mask(bio_curr, bg)
bio_fut <- terra::mask(bio_fut, bg)


#we then need to clip the extent to match an outline of Switzerland
names(bio_curr) <- c('bio_2', 'bio_5', 'bio_14')
names(bio_fut) <- c('bio_2', 'bio_5', 'bio_14')

#plot the current climate variables and how they are projected to change over time
plot(bio_curr)

plot(bio_fut)

#3. Model fitting 
#can you code a binomial GLM with all three predictors fitted as linear and squared terms?

model <- glm( Turdus_torquatus ~ bio_2 + I(bio_2^2) + bio_5 + I(bio_5^2) + bio_14 + I(bio_14^2), family='binomial', data=ouzel_df)

summary(model)

#from the model summary, which variables look to be having an important effect?
#bio5 

#4 testing and critiquing our model
#4.1 partial effects 
#model critique is to understand how each climate variable affects the presence/absence of ring ouzels
#we will plot the partial effects, i.e. the effects of each variable whilst holding the other variables constant
#predict function for each variable. In order to use the predict function we will generate a new data frame in which our focal variable varies but the other variables are held constant.

#bio2
bio_2<-seq(min(ouzel_df$bio_2),max(ouzel_df$bio_2),0.1)

newdata_bio2<-expand.grid(bio_2,mean(ouzel_df$bio_5),mean(ouzel_df$bio_14))
names(newdata_bio2)<-c("bio_2","bio_5","bio_14")
#To use predict we need to generate a new dataset. We will set the other two climate variables at their mean.

response<-predict(model,newdata=newdata_bio2,type="response")
#We've told the predict function to make a prediction on the response scale, ie. in terms of presence/absence

plot(newdata_bio2$bio_2,response,type="l")

#bio14
bio_14<-seq(min(ouzel_df$bio_14),max(ouzel_df$bio_14),0.1)

newdata_bio14<-expand.grid(bio_14,mean(ouzel_df$bio_2),mean(ouzel_df$bio_5))
names(newdata_bio14)<-c("bio_14","bio_2","bio_5")
#To use predict we need to generate a new dataset. We will set the other two climate variables at their mean.

response<-predict(model,newdata=newdata_bio14,type="response")
#We've told the predict function to make a prediction on the response scale, ie. in terms of presence/absence

plot(newdata_bio14$bio_14,response,type="l")

#bio5
bio_5<-seq(min(ouzel_df$bio_5),max(ouzel_df$bio_5),0.1)

newdata_bio5<-expand.grid(bio_5,mean(ouzel_df$bio_2),mean(ouzel_df$bio_14))
names(newdata_bio5)<-c("bio_5","bio_2","bio_14")
#To use predict we need to generate a new dataset. We will set the other two climate variables at their mean.

response<-predict(model,newdata=newdata_bio5,type="response")
#We've told the predict function to make a prediction on the response scale, ie. in terms of presence/absence

plot(newdata_bio5$bio_5,response,type="l")


#4.2 spatial cross validation
#examine how well model can predict presence/absence
#we randomly drop a block from the dataset. (2) We then refit the glm to this reduced “training” dataset. (3) Use the climate data in the ring ouzel dataset to obtain a predicted probability of ring ouzel presence. (4) Consider different probability thresholds to generate a ROC curve. (5) Calculate the area under the curve. And then repeat for each of the blocks.
#In our case there are only 5 spatial blocks, so we will run the cross validation five times and calculate the AUC each time. From these five blocks we can then calculate an average AUC

#first we will exclude spatial block 1.
training1<-subset(ouzel_df,blockCV_tile!=1)
#next we re-run the glm
model1<-glm( Turdus_torquatus ~ bio_2 + I(bio_2^2) + bio_5 + I(bio_5^2) + bio_14 + I(bio_14^2), family='binomial', data=training1)
#we will then subset the data for a testing dataset and we will see how well the glm fitted to the other data does in predicting presences in this testing block
testing<-subset(ouzel_df,blockCV_tile==1)
predicted<-predict(model1,testing,type="response")
#so here we have a prediction as a proportion rather than as a pres/abs.

#The next step is to take different threshold values (i.e probability values at which we count a species as present or absent)
thresholds<-seq(0,1,0.001)
tpr<-c()
fpr<-c()
#We will then use a loop to consider the true positive and false positive rate at each threshold value
for(x in 1:length(thresholds)){
  
  predicted.pres<-as.numeric(predicted>thresholds[x])
  
  present_correct<-length(which(predicted.pres*testing$Turdus_torquatus==1))
  present_incorrect<-length(which(testing$Turdus_torquatus-predicted.pres==1))
  tpr[x]<-present_correct/(present_correct+present_incorrect)
  
  absent_correct<-length(which(predicted.pres+testing$Turdus_torquatus==0))
  absent_incorrect<-length(which(testing$Turdus_torquatus-predicted.pres==-1))
  fpr[x]<-absent_incorrect/(absent_incorrect+absent_correct)
  
}

#When we've run that we can plot the receiver operating characteristic (ROC) curve
plot(fpr,tpr,xlab="false positive rate, 1-sensitivity",ylab="true positive rate, specificity",type="l")
abline(0,1,lty=2)

#Finally to calculate AUC, we can imagine lots of small rectangles. For each one we will calculate its areas and then the area under curve is the sum of those

sortedvals<-cbind(fpr,tpr)[order(fpr),]

AUC<-0
for(x in 2:length(sortedvals[,1])){
  AUC<-AUC+(sortedvals[x,1]-sortedvals[x-1,1])*sortedvals[x-1,2]}

AUC

#quite a high AUC value, suggesting a good model performance when it comes to predicting ring ouzel presence/absence in another region
#want to decide on an optimal probability threshold, as we’ll use this for the mapping stage that comes next
#Here we'll find a value that maximises true positive rate and minimises false positive rate
bestthreshold<-thresholds[which.max(tpr+(1-fpr))]

bestthreshold


