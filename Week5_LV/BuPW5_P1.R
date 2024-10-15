rm(list=ls())
gc()
graphics.off()

library(tidyverse)
library(ggplot2)
library(deSolve)

#the predator-prey lotka-volterra model
#implementation of the base LV model

LG <- function(t,state,parameters){ ##logistic grown function, that takes a set of parameter values, initial conditions and a time sequence
  with(as.list(c(state, parameters)),{ ##"with" is a function that allows us to use the variable names directly - it looks for r, K and P in state and parameters
    
    dP <- r*(1-P/K)*P ##this is our logistic equation governing the rate of change of P
    
    return(list(dP)) ## return the rate of change - it needs to be a list
  }) # end with(as.list ...
}

#modify this function to include two differential equations representing a predator-prey LV model
#the two equations must be included within the same function
#using ode() function, like last time 

install.packages("deSolve")
library(deSolve)

LV <- function(t,state,parameters){
  with(as.list(c(state, parameters)),{
    # rate of change
    dX <- a*X - b*X*Y
    dY <- c*X*Y - d*Y
    
    # return the rate of change
    list(c(dX, dY))
  }) # end with(as.list ...
}


parameters <- c(a=0.05, b=0.025, c=0.025, d=0.2)
state <- c(X=10, Y=10)
times <- seq(0,500,by=0.01)
out <- ode(y=state, times = times, func = LV, parms = parameters)
out.df <- data.frame(out)

#plotting the output 
install.packages("ggplot")
library(ggplot2)

ggplot(data = out.df)+
  geom_line(mapping=aes(x=time,y=X),color="darkolivegreen") +
  geom_line(mapping=aes(x=time,y=Y),color="maroon") +
  geom_hline(yintercept=0,color="darkgrey") +
  geom_vline(xintercept=0,color="darkgrey") +
  labs(x = "Time", y = "P")

#plotting the changes in populations through time in phase space 

ggplot(data = out.df)+
  geom_path(mapping=aes(x=X,y=Y),color="hotpink1") +
  xlim(0,70) +
  ylim(0,40) +
  geom_hline(yintercept=0,color="darkgrey") +
  geom_vline(xintercept=0,color="darkgrey") +
  labs(x = "Prey", y = "Predator")

#when all parameters increased, phase space circle is smaller 
#when all parameters decresed, phase space larger and less circular 

#prey growth rate: exponential vs logistic 
#in the original LV model, prey growth rate is exponential in absence of predator 
#change the equation to include a logistic growth instead 

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
  geom_line(mapping=aes(x=time,y=X),color="darkolivegreen") +
  geom_line(mapping=aes(x=time,y=Y),color="maroon") +
  geom_hline(yintercept=0,color="darkgrey") +
  geom_vline(xintercept=0,color="darkgrey") +
  labs(x = "Time", y = "P")

ggplot(data = out.df)+
  geom_path(mapping=aes(x=X,y=Y),color="hotpink1") +
  xlim(0,70) +
  ylim(0,40) +
  geom_hline(yintercept=0,color="darkgrey") +
  geom_vline(xintercept=0,color="darkgrey") +
  labs(x = "Prey", y = "Predator")

#incorporating functional response 
#at some point predators will not be able to increase their prey intake further 
#the rate at which predators can consume prey is called functional response 
#three types of functional response, including the liner one (type I) that is captured by terms in the original LV equations 
#A controls the slope of the functional response - high value indicates poor pred hunting efficiency

#A <- 0.005
A <- 0
x <- seq(0,30,0.1)
y <- x/(1+A*x)
ggplot()+
  geom_line(mapping=aes(x=x,y=y),color="darkslategrey") +
  geom_hline(yintercept=0,color="darkgrey") +
  geom_vline(xintercept=0,color="darkgrey") +
  labs(x = "Prey population", y = "Prey consumed")

#include a type II functional response in your implementation of the Lotka-Volterra model
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
  geom_line(mapping=aes(x=time,y=X),color="darkolivegreen") +
  geom_line(mapping=aes(x=time,y=Y),color="maroon") +
  geom_hline(yintercept=0,color="darkgrey") +
  geom_vline(xintercept=0,color="darkgrey") +
  labs(x = "Time", y = "P")

ggplot(data = out.df)+
  geom_path(mapping=aes(x=X,y=Y),color="hotpink1") +
  xlim(0,270) +
  ylim(0,150) +
  geom_hline(yintercept=0,color="darkgrey") +
  geom_vline(xintercept=0,color="darkgrey") +
  labs(x = "Prey", y = "Predator")

#combine both the logistic growth and the functional response in your Lotka-Volterra model

x <- seq(0,30,0.1)
A <- 0.1
y <- x/(1+A*x)
ggplot()+
  geom_line(mapping=aes(x=x,y=y),color="darkslategrey") +
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
  geom_line(mapping=aes(x=time,y=X),color="darkolivegreen") +
  geom_line(mapping=aes(x=time,y=Y),color="maroon") +
  geom_hline(yintercept=0,color="darkgrey") +
  geom_vline(xintercept=0,color="darkgrey") +
  labs(x = "Time", y = "P")


#three-species competition Lotka-Volterra model: limiting similarity

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

ttimes <- seq(0,1000,by=0.01)

out <- ode(y=state, times = times, func = LS, parms = parameters)

out.df <- data.frame(out)

ggplot(data = out.df)+
  geom_line(mapping=aes(x=time,y=X1),color="darkolivegreen") +
  geom_line(mapping=aes(x=time,y=X2),color="maroon") +
  geom_hline(yintercept=0,color="darkgrey") +
  geom_vline(xintercept=0,color="darkgrey") +
  labs(x = "Time", y = "P")


ggplot(data = out.df)+
  geom_path(mapping=aes(x=X1,y=X2),color="hotpink1") +
  geom_hline(yintercept=0,color="darkgrey") +
  geom_vline(xintercept=0,color="darkgrey") +
  labs(x = "Species 1", y = "Species 2")

#adding species 3
#aij parameters will be determined by the relative niched of the specie
alpha.func <- function(mu1,sig1,mu2,sig2,K1,K2,start,end){ ##this is the function to compute the alpha coefficients from the mean and standard deviations of the Gaussian niches of the species and the start and end values of the environment
  niche1 <- K1*dnorm(seq(start,end,length.out=100),mean=mu1,sd=sig1) ##dnorm() generates the values of the Gaussian. Check ?dnorm
  niche2 <- K2*dnorm(seq(start,end,length.out=100),mean=mu2,sd=sig2)
  a <- sum(niche1*niche2)/sum(niche1*niche1) ##because we have discrete values, we use a sum to approximate the integral
  return(a)
}

#trying different parameter values
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

#visualising the niches 

#setup and solve the system of differential equations
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

#plot the populations
ggplot(data = out.df)+
  geom_line(mapping=aes(x=time,y=X1),color="darkolivegreen") +
  geom_line(mapping=aes(x=time,y=X2),color="maroon") +
  geom_line(mapping=aes(x=time,y=X3),color="darkslategrey") +
  geom_hline(yintercept=0,color="darkgrey") +
  geom_vline(xintercept=0,color="darkgrey") +
  labs(x = "Time", y = "P")

