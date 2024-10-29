library(tidyverse)
library(vegan)

install.packages("betapart")
library(betapart)

dat.all <- read.csv("Species_list_full_2905.csv",sep=";")
dat.all <- dat.all[-which(is.na(dat.all$species)),]
dat.all <- dat.all[-which(is.na(dat.all$island)),]

which(is.na(dat.all$species)) ##which species have NA’s
which(is.na(dat.all$island)) ##which islands have NA’s

#creating 6 data frames 
dat.soc <- dat.all[which(dat.all$islandgroup=="Society"),]
dat.haw <- dat.all[which(dat.all$islandgroup=="Hawaiian"),]
dat.sam <- dat.all[which(dat.all$islandgroup=="Samoa"),]
dat.mar <- dat.all[which(dat.all$islandgroup=="Marquesas"),]
dat.fij <- dat.all[which(dat.all$islandgroup=="Fiji"),]
dat.comb <- rbind(dat.soc,dat.haw,dat.sam,dat.mar,dat.fij)

#computing community patterns
#To compute community patterns, we will need site-by-species data frames. We need to create these data frames as follows, one for each tidy data frame:
#•	create temporary data frames with 3 columns: one for islands, one for species, and an extra one filled with 1’s to represent presences
#•	use the function pivot_wider()
#•	replace NA’s by 0’s
#•	rename the rows with the island names stored in the first column
#•	delete the first column

dat.soc.red <- dat.soc[,c("species","island")] 
dat.soc.red$presence <- 1 
dat.soc.pa <- dat.soc.red %>% 
  pivot_wider(names_from=species,values_from=c(presence))
list0 <- as.list(rep(0,ncol(dat.soc.pa)))
names(list0) <- names(dat.soc.pa)
dat.soc.pa <- as.data.frame(dat.soc.pa %>% replace_na(list0))
row.names(dat.soc.pa) <- dat.soc.pa$island
dat.soc.pa <- dat.soc.pa[,-1]

dat.haw.red <- dat.haw[,c("species","island")] 
dat.haw.red$presence <- 1 
dat.haw.pa <- dat.haw.red %>% 
  pivot_wider(names_from=species,values_from=c(presence))
list0 <- as.list(rep(0,ncol(dat.haw.pa)))
names(list0) <- names(dat.haw.pa)
dat.haw.pa <- as.data.frame(dat.haw.pa %>% replace_na(list0))
row.names(dat.haw.pa) <- dat.haw.pa$island
dat.haw.pa <- dat.haw.pa[,-1]

dat.sam.red <- dat.sam[,c("species","island")] 
dat.sam.red$presence <- 1 
dat.sam.pa <- dat.sam.red %>% 
  pivot_wider(names_from=species,values_from=c(presence))
list0 <- as.list(rep(0,ncol(dat.sam.pa)))
names(list0) <- names(dat.sam.pa)
dat.sam.pa <- as.data.frame(dat.sam.pa %>% replace_na(list0))
row.names(dat.sam.pa) <- dat.sam.pa$island
dat.sam.pa <- dat.sam.pa[,-1]

dat.mar.red <- dat.mar[,c("species","island")] 
dat.mar.red$presence <- 1 
dat.mar.pa <- dat.mar.red %>% 
  pivot_wider(names_from=species,values_from=c(presence))
list0 <- as.list(rep(0,ncol(dat.mar.pa)))
names(list0) <- names(dat.mar.pa)
dat.mar.pa <- as.data.frame(dat.mar.pa %>% replace_na(list0))
row.names(dat.mar.pa) <- dat.mar.pa$island
dat.mar.pa <- dat.mar.pa[,-1]

dat.fij.red <- dat.fij[,c("species","island")] 
dat.fij.red$presence <- 1 
dat.fij.pa <- dat.fij.red %>% 
  pivot_wider(names_from=species,values_from=c(presence))
list0 <- as.list(rep(0,ncol(dat.fij.pa)))
names(list0) <- names(dat.fij.pa)
dat.fij.pa <- as.data.frame(dat.fij.pa %>% replace_na(list0))
row.names(dat.fij.pa) <- dat.fij.pa$island
dat.fij.pa <- dat.fij.pa[,-1]

dat.comb.red <- dat.comb[,c("species","island")] 
dat.comb.red$presence <- 1 
dat.comb.pa <- dat.comb.red %>% 
  pivot_wider(names_from=species,values_from=c(presence))
list0 <- as.list(rep(0,ncol(dat.comb.pa)))
names(list0) <- names(dat.comb.pa)
dat.comb.pa <- as.data.frame(dat.comb.pa %>% replace_na(list0))
row.names(dat.comb.pa) <- dat.comb.pa$island
dat.comb.pa <- dat.comb.pa[,-1]

#gamma and alpha diversity 
#gamma is just ncol, because it is a measure of the total species diversity 
ncol(dat.soc.pa) #684
ncol(dat.haw.pa) #1770 <-
ncol(dat.sam.pa) #412
ncol(dat.mar.pa) #474
ncol(dat.fij.pa) #533
ncol(dat.comb.pa) #2069

#alpha diversity
mean(rowSums(dat.soc.pa)) #149.8571
mean(rowSums(dat.haw.pa)) #416.4615 <- 
mean(rowSums(dat.sam.pa)) #109.3333
mean(rowSums(dat.mar.pa)) #165.6
mean(rowSums(dat.fij.pa)) #30.7551
mean(rowSums(dat.comb.pa)) #122.3163

#species accumulation curves 
SAC.soc <- specaccum(dat.soc.pa)
SAC.haw <- specaccum(dat.haw.pa)
SAC.sam <- specaccum(dat.sam.pa)
SAC.mar <- specaccum(dat.mar.pa)
SAC.fij <- specaccum(dat.fij.pa)
SAC.comb <- specaccum(dat.comb.pa)

Estim.soc <- poolaccum(dat.soc.pa)
Estim.haw <- poolaccum(dat.haw.pa)
Estim.sam <- poolaccum(dat.sam.pa)
Estim.mar <- poolaccum(dat.mar.pa)
Estim.fij <- poolaccum(dat.fij.pa)
Estim.comb <- poolaccum(dat.comb.pa)

SAC.soc 
SAC.haw 
SAC.sam 
SAC.mar 
SAC.fij 
SAC.comb 

par(mfrow=c(2,3)) #allows two rows of 3 graphs 
plot(SAC.soc$richness,pch=1,lty=1,lwd=2,type="b",col="darkolivegreen",ylim=c(0,max(rowMeans(Estim.soc$chao))),ylab="Richness",main="Society")
points(3:nrow(dat.soc.pa),rowMeans(Estim.soc$chao),pch=2,lty=2,lwd=2,type="b",col="maroon")
plot(SAC.haw$richness,pch=1,lty=1,lwd=2,type="b",col="darkolivegreen",ylim=c(0,max(rowMeans(Estim.haw$chao))),ylab="Richness",main="Hawai'i")
points(3:nrow(dat.haw.pa),rowMeans(Estim.haw$chao),pch=2,lty=2,lwd=2,type="b",col="maroon")
plot(SAC.sam$richness,pch=1,lty=1,lwd=2,type="b",col="darkolivegreen",ylim=c(0,max(rowMeans(Estim.sam$chao))),ylab="Richness",main="Samoa")
points(3:nrow(dat.sam.pa),rowMeans(Estim.sam$chao),pch=2,lty=2,lwd=2,type="b",col="maroon")
plot(SAC.mar$richness,pch=1,lty=1,lwd=2,type="b",col="darkolivegreen",ylim=c(0,max(rowMeans(Estim.mar$chao))),ylab="Richness",main="Marquesas")
points(3:nrow(dat.mar.pa),rowMeans(Estim.mar$chao),pch=2,lty=2,lwd=2,type="b",col="maroon")
plot(SAC.fij$richness,pch=1,lty=1,lwd=2,type="b",col="darkolivegreen",ylim=c(0,max(rowMeans(Estim.fij$chao))),ylab="Richness",main="Fiji")
points(3:nrow(dat.fij.pa),rowMeans(Estim.fij$chao),pch=2,lty=2,lwd=2,type="b",col="maroon")
plot(SAC.comb$richness,pch=1,lty=1,lwd=2,type="b",col="darkolivegreen",ylim=c(0,max(rowMeans(Estim.comb$chao))),ylab="Richness",main="All data")
points(3:nrow(dat.comb.pa),rowMeans(Estim.comb$chao),pch=2,lty=2,lwd=2,type="b",col="maroon")

#beta diversity 
beta.soc <- beta.pair(dat.soc.pa)
beta.haw <- beta.pair(dat.haw.pa)
beta.sam <- beta.pair(dat.sam.pa)
beta.mar <- beta.pair(dat.mar.pa)
beta.fij <- beta.pair(dat.fij.pa)
beta.comb <- beta.pair(dat.comb.pa)

mean(beta.soc$beta.sim) #0.2779589
mean(beta.haw$beta.sim) #0.176857
mean(beta.sam$beta.sim) #0.4066672
mean(beta.mar$beta.sim) #0.1442509 <- most similar (accounting for spatial turnover)
mean(beta.fij$beta.sim) #0.7408851
mean(beta.comb$beta.sim) #0.6095149

mean(beta.soc$beta.sor) #0.7272671
mean(beta.haw$beta.sor) #7896844
mean(beta.sam$beta.sor) #0.7302065
mean(beta.mar$beta.sor) #0.5947258 <- most similar (total dissimilarity)
mean(beta.fij$beta.sor) #0.9384918
mean(beta.comb$beta.sor) #0.9083158

#plotting the islands based on dissimilarity 

coord.comb.sim <- data.frame(cmdscale(beta.comb$beta.sim))
coord.comb.sim$col <- c(rep("blue",nrow(dat.soc.pa)),rep("red",nrow(dat.haw.pa)),rep("darkgreen",nrow(dat.sam.pa)),rep("orange",nrow(dat.mar.pa)),rep("purple",nrow(dat.fij.pa)))
coord.comb.sim$pch <- c(rep(0,nrow(dat.soc.pa)),rep(1,nrow(dat.haw.pa)),rep(2,nrow(dat.sam.pa)),rep(3,nrow(dat.mar.pa)),rep(4,nrow(dat.fij.pa)))

coord.comb.sor <- data.frame(cmdscale(beta.comb$beta.sor))
coord.comb.sor$col <- c(rep("blue",nrow(dat.soc.pa)),rep("red",nrow(dat.haw.pa)),rep("darkgreen",nrow(dat.sam.pa)),rep("orange",nrow(dat.mar.pa)),rep("purple",nrow(dat.fij.pa)))
coord.comb.sor$pch <- c(rep(0,nrow(dat.soc.pa)),rep(1,nrow(dat.haw.pa)),rep(2,nrow(dat.sam.pa)),rep(3,nrow(dat.mar.pa)),rep(4,nrow(dat.fij.pa)))

par(mfrow=c(1,3))
plot(coord.comb.sim$X1,coord.comb.sim$X2,col=coord.comb.sim$col,pch=coord.comb.sim$pch,lwd=2,main="Simpson")
plot(coord.comb.sor$X1,coord.comb.sor$X2,col=coord.comb.sor$col,pch=coord.comb.sor$pch,lwd=2,main="Sorensen")
plot.new()
legend(x="topleft",legend=c("Society","Hawai’i","Samoa","Marquesas","Fiji"),col = c("blue","red","darkgreen","orange","purple"),pch=0:4,lwd=2,bty="n",lty=0,cex=2)

#occupation frequency distributions
freq.soc <- colSums(dat.soc.pa)/nrow(dat.soc.pa)
freq.haw <- colSums(dat.haw.pa)/nrow(dat.haw.pa)
freq.sam <- colSums(dat.sam.pa)/nrow(dat.sam.pa)
freq.mar <- colSums(dat.mar.pa)/nrow(dat.mar.pa)
freq.fij <- colSums(dat.fij.pa)/nrow(dat.fij.pa)
freq.comb <- colSums(dat.comb.pa)/nrow(dat.comb.pa)

par(mfrow=c(2,3))
hist(freq.soc,main="Society",xlab="Occupancy",breaks = seq(0,1,0.1))
hist(freq.haw,main="Hawa'i",xlab="Occupancy",breaks = seq(0,1,0.1))
hist(freq.sam,main="Samoa",xlab="Occupancy",breaks = seq(0,1,0.1))
hist(freq.mar,main="Marquesas",xlab="Occupancy",breaks = seq(0,1,0.1))
hist(freq.fij,main="Fiji",xlab="Occupancy",breaks = seq(0,1,0.1))
hist(freq.comb,main="Combined",xlab="Occupancy",breaks = seq(0,1,0.1))

