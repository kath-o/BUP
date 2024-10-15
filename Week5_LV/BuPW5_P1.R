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


parameters <- c(a=0.1, b=0.02, c=0.02, d=0.4)
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


