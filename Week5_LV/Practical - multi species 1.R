rm(list=ls())
gc()
graphics.off()

library(tidyverse)
library(ggplot2)
library(deSolve)


#################
##Predator-Prey##
#################

##original LV

LV <- function(t,state,parameters){
  with(as.list(c(state, parameters)),{
    # rate of change
    dX <- a*X - b*X*Y
    dY <- c*X*Y - d*Y
    
    # return the rate of change
    list(c(dX, dY))
  }) # end with(as.list ...
}

parameters <- c(a=0.1, b=0.02, c=0.02, d=0.4)
state <- c(X=10, Y=10)
times <- seq(0,500,by=0.01)
out <- ode(y=state, times = times, func = LV, parms = parameters)
out.df <- data.frame(out)

ggplot(data = out.df)+
  geom_line(mapping=aes(x=time,y=X),color="blue") +
  geom_line(mapping=aes(x=time,y=Y),color="red") +
  geom_hline(yintercept=0,color="darkgrey") +
  geom_vline(xintercept=0,color="darkgrey") +
  labs(x = "Time", y = "P")


ggplot(data = out.df)+
  geom_path(mapping=aes(x=X,y=Y),color="red") +
  xlim(0,70) +
  ylim(0,40) +
  geom_hline(yintercept=0,color="darkgrey") +
  geom_vline(xintercept=0,color="darkgrey") +
  labs(x = "Prey", y = "Predator")



##LV with logistic growth


LV.lg <- function(t,state,parameters){
  with(as.list(c(state, parameters)),{
    # rate of change
    dX <- a*X*(1-X/K) - b*X*Y
    dY <- c*X*Y - d*Y
    
    # return the rate of change
    list(c(dX, dY))
  }) # end with(as.list ...
}

parameters <- c(a=0.1, b=0.02, c=0.02, d=0.4, K=30)
state <- c(X=10, Y=10)
times <- seq(0,500,by=0.01)
out <- ode(y=state, times = times, func = LV.lg, parms = parameters)
out.df <- data.frame(out)

ggplot(data = out.df)+
  geom_line(mapping=aes(x=time,y=X),color="blue") +
  geom_line(mapping=aes(x=time,y=Y),color="red") +
  geom_hline(yintercept=0,color="darkgrey") +
  geom_vline(xintercept=0,color="darkgrey") +
  labs(x = "Time", y = "P")


ggplot(data = out.df)+
  geom_path(mapping=aes(x=X,y=Y),color="red") +
  xlim(0,70) +
  ylim(0,40) +
  geom_hline(yintercept=0,color="darkgrey") +
  geom_vline(xintercept=0,color="darkgrey") +
  labs(x = "Prey", y = "Predator")



##LV with functional response

#A <- 0.005
A <- 0
x <- seq(0,30,0.1)
y <- x/(1+A*x)
ggplot()+
  geom_line(mapping=aes(x=x,y=y),color="blue") +
  geom_hline(yintercept=0,color="darkgrey") +
  geom_vline(xintercept=0,color="darkgrey") +
  labs(x = "Prey population", y = "Prey consumed")


LV.fr <- function(t,state,parameters){
  with(as.list(c(state, parameters)),{
    # rate of change
    dX <- a*X - b*X*Y/(1+A*X)
    dY <- c*X*Y/(1+A*X) - d*Y
    
    # return the rate of change
    list(c(dX, dY))
  }) # end with(as.list ...
}

parameters <- c(a=0.1, b=0.02, c=0.02, d=0.4, A=0.02)
state <- c(X=10, Y=10)
times <- seq(0,500,by=0.01)
out <- ode(y=state, times = times, func = LV.fr, parms = parameters)
out.df <- data.frame(out)


ggplot(data = out.df)+
  geom_line(mapping=aes(x=time,y=X),color="blue") +
  geom_line(mapping=aes(x=time,y=Y),color="red") +
  geom_hline(yintercept=0,color="darkgrey") +
  geom_vline(xintercept=0,color="darkgrey") +
  labs(x = "Time", y = "P")


ggplot(data = out.df)+
  geom_path(mapping=aes(x=X,y=Y),color="red") +
  xlim(0,270) +
  ylim(0,150) +
  geom_hline(yintercept=0,color="darkgrey") +
  geom_vline(xintercept=0,color="darkgrey") +
  labs(x = "Prey", y = "Predator")



##LV with logistic growth and functional response

x <- seq(0,30,0.1)
A <- 0.1
y <- x/(1+A*x)
ggplot()+
  geom_line(mapping=aes(x=x,y=y),color="blue") +
  geom_hline(yintercept=0,color="darkgrey") +
  geom_vline(xintercept=0,color="darkgrey") +
  labs(x = "Prey population", y = "Prey consumed")


LV.lg.fr <- function(t,state,parameters){
  with(as.list(c(state, parameters)),{
    # rate of change
    dX <- a*X*(1-X/K) - b*X*Y/(1+A*X)
    dY <- c*X*Y/(1+A*X) - d*Y
    
    # return the rate of change
    list(c(dX, dY))
  }) # end with(as.list ...
}

parameters <- c(a=0.1, b=0.02, c=0.02, d=0.4, K=30, A=0.01)
state <- c(X=10, Y=10)
times <- seq(0,500,by=0.01)
out <- ode(y=state, times = times, func = LV.lg.fr, parms = parameters)
out.df <- data.frame(out)


ggplot(data = out.df)+
  geom_line(mapping=aes(x=time,y=X),color="blue") +
  geom_line(mapping=aes(x=time,y=Y),color="red") +
  geom_hline(yintercept=0,color="darkgrey") +
  geom_vline(xintercept=0,color="darkgrey") +
  labs(x = "Time", y = "P")


ggplot(data = out.df)+
  geom_path(mapping=aes(x=X,y=Y),color="red") +
  xlim(0,70) +
  ylim(0,40) +
  geom_hline(yintercept=0,color="darkgrey") +
  geom_vline(xintercept=0,color="darkgrey") +
  labs(x = "Prey", y = "Predator")




#################
##Competition##
#################


##2 species

parameters <- c(a12=1, a21=0.9, r=0.3, K = 100)
state <- c(X1=50, X2=10)

LS <- function(t,state,parameters){
  with(as.list(c(state, parameters)),{
    # rate of change
    dX1 <- r*X1*(1-(X1+a12*X2)/K)
    dX2 <- r*X2*(1-(X2+a21*X1)/K)
    
    # return the rate of change
    list(c(dX1, dX2))
  }) # end with(as.list ...
}

times <- seq(0,1000,by=0.01)

out <- ode(y=state, times = times, func = LS, parms = parameters)

out.df <- data.frame(out)

ggplot(data = out.df)+
  geom_line(mapping=aes(x=time,y=X1),color="blue") +
  geom_line(mapping=aes(x=time,y=X2),color="red") +
  geom_hline(yintercept=0,color="darkgrey") +
  geom_vline(xintercept=0,color="darkgrey") +
  labs(x = "Time", y = "P")


ggplot(data = out.df)+
  geom_path(mapping=aes(x=X1,y=X2),color="red") +
  geom_hline(yintercept=0,color="darkgrey") +
  geom_vline(xintercept=0,color="darkgrey") +
  labs(x = "Species 1", y = "Species 2")

##3 species

alpha.func <- function(mu1,sig1,mu2,sig2,K1,K2,start,end){ ##this is the function to compute the alpha coefficients from the mean and standard deviations of the Gaussian niches of the species and the start and end values of the environment
  niche1 <- K1*dnorm(seq(start,end,length.out=100),mean=mu1,sd=sig1) ##dnorm() generates the values of the Gaussian. Check ?dnorm
  niche2 <- K2*dnorm(seq(start,end,length.out=100),mean=mu2,sd=sig2)
  a <- sum(niche1*niche2)/sum(niche1*niche1) ##because we have discrete values, we use a sum to approximate the integral
  return(a)
}

##Let's try different parameter values
D <- 10 ##distance between the niche optima
mu1 <- 5 ##niche optima of species 1
mu2 <- mu1+D ##niche optima of species 2
mu3 <- mu1+2*D ##niche optima of species 3
sig1 <- sig2 <- sig3 <- 10 ##all species niches have the same standard deviation for simplicity
start <- 0
end <- 30
K1 <- 200 ##carrying capacity species 1 and 3
K2 <- 250 ##carrying capacity species 2
a12 <- alpha.func(mu1,sig1,mu2,sig2,K1,K2,start,end)
a13 <- alpha.func(mu1,sig1,mu3,sig3,K1,K1,start,end)
a21 <- alpha.func(mu2,sig2,mu1,sig1,K2,K1,start,end)
a23 <- alpha.func(mu2,sig2,mu3,sig3,K2,K1,start,end)
a31 <- alpha.func(mu3,sig3,mu1,sig1,K1,K1,start,end)
a32 <- alpha.func(mu3,sig3,mu2,sig2,K1,K2,start,end)


##visualise the niches
resource <- seq(start,end,length.out=100)
niche1 <- dnorm(resource,mean=mu1,sd=sig1)*K1
niche2 <- dnorm(resource,mean=mu2,sd=sig2)*K2
niche3 <- dnorm(resource,mean=mu3,sd=sig3)*K1
ggplot()+
  geom_line(mapping=aes(x=resource,y=niche1),color="blue")+
  geom_line(mapping=aes(x=resource,y=niche2),color="red")+
  geom_line(mapping=aes(x=resource,y=niche3),color="darkgreen")


##setup and solve the system of differential equations
parameters <- c(a12=a12, a13=a13, a21=a21, a23=a23, a31=a31, a32=a32, r=0.3, K1 = K1, K2 = K2)
state <- c(X1=10, X2=10, X3=10)

LS2 <- function(t,state,parameters){
  with(as.list(c(state, parameters)),{
    # rate of change
    dX1 <- r*X1*(1-(X1+a12*X2+a13*X3)/K1)
    dX2 <- r*X2*(1-(X2+a21*X1+a23*X3)/K2)
    dX3 <- r*X3*(1-(X3+a31*X1+a32*X2)/K1)
    
    # return the rate of change
    list(c(dX1, dX2, dX3))
  }) # end with(as.list ...
}

times <- seq(0,200,by=0.01)
out <- ode(y=state, times = times, func = LS2, parms = parameters)
out.df <- data.frame(out)

##plot the populations
ggplot(data = out.df)+
  geom_line(mapping=aes(x=time,y=X1),color="blue") +
  geom_line(mapping=aes(x=time,y=X2),color="red") +
  geom_line(mapping=aes(x=time,y=X3),color="darkgreen") +
  geom_hline(yintercept=0,color="darkgrey") +
  geom_vline(xintercept=0,color="darkgrey") +
  labs(x = "Time", y = "P")




