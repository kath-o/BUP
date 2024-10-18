#discrete mechanistic model for 6 populations
#lichen, plants, moose forage, caribou, moose, wolves

#setting parameters for each equation
sP <- 0.5
Kp <- 240
lP <- 2.68275
hp <- 300
sL <- 0.07
Kl <- 870
lL <- 2.68275
hl <- 400
sH <- 0.5
Kh <- 970
lH <- 3.1682
hh <- 1000
fC <- 1
eC <- 1.85
mC <- 0
fM <- 1.5
eM <- 0.6
b <- 0.8
g <- 0.38
dC <- 460
dM <- 46
mM <- 0
mW <- 0

nsteps <- 200
pop.df.1 <- data.frame(time=1:nsteps,P=numeric(nsteps),L=numeric(nsteps),H=numeric(nsteps),C=numeric(nsteps),M=numeric(nsteps),W=numeric(nsteps))

pop.df.1 <- within(pop.df.1,{
  P[1] <- 240
  L[1] <- 870
  H[1] <- 970
  C[1] <- 7
  M[1] <- 25
  W[1] <- 8
})


for(t in 2:nsteps){
  pop.df.1 <- within(pop.df.1,{
    ##Primary producers
    P.birth <- sP*P[t-1]*(1-P[t-1]/Kp)
    P.death <- lP*P[t-1]*C[t-1]/(hp+P[t-1])
    P[t] <- max(0,P[t-1] + P.birth - P.death) ##plants consumed by Caribou - the max is used to avoid negative values
    
    L.birth <- sL*L[t-1]*(1-L[t-1]/Kl)
    L.death <- lL*L[t-1]*C[t-1]/(hl+L[t-1])
    L[t] <- max(0,L[t-1] + L.birth - L.death) ##lichen consumed by Caribou
    
    h.birth <- sH*H[t-1]*(1-H[t-1]/Kh)
    h.death <- lH*H[t-1]*M[t-1]/(hH+H[t-1])
    H[t] <- max(0,H[t-1] + h.birth - h.death) ##plants consumed by Moose
    
    ##First trophic level
    C.growth <- C[t-1]*fC*(P[t-1]/(P[t-1]+hp))*(L[t-1]/(L[t-1]+hl))
    C.death <- C[t-1]*(1-P[t-1]/(P[t-1]+hp))*(1-L[t-1]/(L[t-1]+hl))
    C.pred <-  W[t-1]*eC*C[t-1]/(C[t-1]+dC) 
    C[t] <- max(0,C[t-1] + C.growth - C.death - C.pred)
    
    M.growth <- M[t-1]*fM*H[t-1]/(H[t-1]+hH)
    M.death <- M[t-1]*(1-H[t-1]/(H[t-1]+hH))
    M.pred <-  W[t-1]*eM*M[t-1]/(M[t-1]+dM)
    M[t] <- max(0,M[t-1] + M.growth - M.death - M.pred)
    
    
    ##Predator
    W.growth <- (b*(C[t-1]/(dC+C[t-1])) + b*(M[t-1]/(dM+M[t-1])))*W[t-1]
    W.death <- g*W[t-1]
    W[t] <- max(0,W[t-1] + W.growth - W.death)
  })
}

colors <- c("Plants"="brown","Lichen"="orange","Shrubs"="green","Caribou"="purple","Moose"="blue","Wolf"="red")
ggplot(data = pop.df.1)+
  geom_line(mapping=aes(x=time,y=P,color="Plants")) +
  geom_line(mapping=aes(x=time,y=L,color="Lichen")) +
  geom_line(mapping=aes(x=time,y=H,color="Shrubs")) +
  geom_line(mapping=aes(x=time,y=C,color="Caribou")) +
  geom_line(mapping=aes(x=time,y=M,color="Moose")) +
  geom_line(mapping=aes(x=time,y=W,color="Wolf")) +
  geom_hline(yintercept=0,color="darkgrey") +
  geom_vline(xintercept=0,color="darkgrey") +
  labs(x = "Time", y = "Densities", color="Legend",title="Original model")+
  scale_color_manual(values = colors)

#limit on the number of wolves

sP <- 0.5
Kp <- 240
lP <- 2.68275
hp <- 300
sL <- 0.07
Kl <- 870
lL <- 2.68275
hl <- 400
sH <- 0.5
Kh <- 970
lH <- 3.1682
hH <- 1000
fC <- 1
eC <- 1.85
mC <- 0
fM <- 1.5
eM <- 0.6
b <- 0.8
g <- 0.38
dC <- 460
dM <- 46


nsteps <- 200
pop.df.2 <- data.frame(time=1:nsteps,P=numeric(nsteps),L=numeric(nsteps),H=numeric(nsteps),C=numeric(nsteps),M=numeric(nsteps),W=numeric(nsteps))

pop.df.2 <- within(pop.df.2,{
  P[1] <- 240
  L[1] <- 870
  H[1] <- 970
  C[1] <- 7
  M[1] <- 25
  W[1] <- 8
})


for(t in 2:nsteps){
  pop.df.2 <- within(pop.df.2,{
    ##Primary producers
    P.birth <- sP*P[t-1]*(1-P[t-1]/Kp)
    P.death <- lP*P[t-1]*C[t-1]/(hp+P[t-1])
    P[t] <- max(0,P[t-1] + P.birth - P.death) ##plants consumed by Caribou - the max is used to avoid negative values
    
    L.birth <- sL*L[t-1]*(1-L[t-1]/Kl)
    L.death <- lL*L[t-1]*C[t-1]/(hl+L[t-1])
    L[t] <- max(0,L[t-1] + L.birth - L.death) ##lichen consumed by Caribou
    
    P.birth <- sH*H[t-1]*(1-H[t-1]/Kh)
    P.death <- lH*H[t-1]*M[t-1]/(hH+H[t-1])
    H[t] <- max(0,H[t-1] + P.birth - P.death) ##plants consumed by Moose
    
    ##First trophic level
    C.growth <- C[t-1]*fC*(P[t-1]/(P[t-1]+hp))*(L[t-1]/(L[t-1]+hl))
    C.death <- C[t-1]*(1-P[t-1]/(P[t-1]+hp))*(1-L[t-1]/(L[t-1]+hl))
    C.pred <-  W[t-1]*eC*C[t-1]/(C[t-1]+dC) 
    C[t] <- max(0,C[t-1] + C.growth - C.death - C.pred)
    
    M.growth <- M[t-1]*fM*H[t-1]/(H[t-1]+hH)
    M.death <- M[t-1]*(1-H[t-1]/(H[t-1]+hH))
    M.pred <-  W[t-1]*eM*M[t-1]/(M[t-1]+dM)
    M[t] <- max(0,M[t-1] + M.growth - M.death - M.pred)
    
    
    ##Predator
    W.growth <- (b*(C[t-1]/(dC+C[t-1])) + b*(M[t-1]/(dM+M[t-1])))*W[t-1]
    W.death <- g*W[t-1]
    W[t] <- min(max(0,W[t-1] + W.growth - W.death),10)
  })
}


colors <- c("Plants"="brown","Lichen"="orange","Shrubs"="green","Caribou"="purple","Moose"="blue","Wolf"="red")
ggplot(data = pop.df.2)+
  geom_line(mapping=aes(x=time,y=P,color="Plants")) +
  geom_line(mapping=aes(x=time,y=L,color="Lichen")) +
  geom_line(mapping=aes(x=time,y=H,color="Shrubs")) +
  geom_line(mapping=aes(x=time,y=C,color="Caribou")) +
  geom_line(mapping=aes(x=time,y=M,color="Moose")) +
  geom_line(mapping=aes(x=time,y=W,color="Wolf")) +
  geom_hline(yintercept=0,color="darkgrey") +
  geom_vline(xintercept=0,color="darkgrey") +
  labs(x = "Time", y = "Densities", color="Legend",title="Wolf cap")+
  scale_color_manual(values = colors)

#model with limit on moose population 
sP <- 0.5
Kp <- 240
lP <- 2.68275
hp <- 300
sL <- 0.07
Kl <- 870
lL <- 2.68275
hl <- 400
sH <- 0.5
Kh <- 970
lH <- 3.1682
hH <- 1000
fC <- 1
eC <- 1.85
mC <- 0
fM <- 1.5
eM <- 0.6
b <- 0.8
g <- 0.38
dC <- 460
dM <- 46


nsteps <- 200
pop.df.3 <- data.frame(time=1:nsteps,P=numeric(nsteps),L=numeric(nsteps),H=numeric(nsteps),C=numeric(nsteps),M=numeric(nsteps),W=numeric(nsteps))

pop.df.3 <- within(pop.df.3,{
  P[1] <- 240
  L[1] <- 870
  H[1] <- 970
  C[1] <- 7
  M[1] <- 25
  W[1] <- 8
})


for(t in 2:nsteps){
  pop.df.3 <- within(pop.df.3,{
    ##Primary producers
    P.birth <- sP*P[t-1]*(1-P[t-1]/Kp)
    P.death <- lP*P[t-1]*C[t-1]/(hp+P[t-1])
    P[t] <- max(0,P[t-1] + P.birth - P.death) ##plants consumed by Caribou - the max is used to avoid negative values
    
    L.birth <- sL*L[t-1]*(1-L[t-1]/Kl)
    L.death <- lL*L[t-1]*C[t-1]/(hl+L[t-1])
    L[t] <- max(0,L[t-1] + L.birth - L.death) ##lichen consumed by Caribou
    
    P.birth <- sH*H[t-1]*(1-H[t-1]/Kh)
    P.death <- lH*H[t-1]*M[t-1]/(hH+H[t-1])
    H[t] <- max(0,H[t-1] + P.birth - P.death) ##plants consumed by Moose
    
    ##First trophic level
    C.growth <- C[t-1]*fC*(P[t-1]/(P[t-1]+hp))*(L[t-1]/(L[t-1]+hl))
    C.death <- C[t-1]*(1-P[t-1]/(P[t-1]+hp))*(1-L[t-1]/(L[t-1]+hl))
    C.pred <-  W[t-1]*eC*C[t-1]/(C[t-1]+dC) 
    C[t] <- max(0,C[t-1] + C.growth - C.death - C.pred)
    
    M.growth <- M[t-1]*fM*H[t-1]/(H[t-1]+hH)
    M.death <- M[t-1]*(1-H[t-1]/(H[t-1]+hH))
    M.pred <-  W[t-1]*eM*M[t-1]/(M[t-1]+dM)
    M[t] <- min(max(0,M[t-1] + M.growth - M.death - M.pred),30)
    
    
    ##Predator
    W.growth <- (b*(C[t-1]/(dC+C[t-1])) + b*(M[t-1]/(dM+M[t-1])))*W[t-1]
    W.death <- g*W[t-1]
    W[t] <- max(0,W[t-1] + W.growth - W.death)
  })
}


colors <- c("Plants"="brown","Lichen"="orange","Shrubs"="green","Caribou"="purple","Moose"="blue","Wolf"="red")
ggplot(data = pop.df.3)+
  geom_line(mapping=aes(x=time,y=P,color="Plants")) +
  geom_line(mapping=aes(x=time,y=L,color="Lichen")) +
  geom_line(mapping=aes(x=time,y=H,color="Shrubs")) +
  geom_line(mapping=aes(x=time,y=C,color="Caribou")) +
  geom_line(mapping=aes(x=time,y=M,color="Moose")) +
  geom_line(mapping=aes(x=time,y=W,color="Wolf")) +
  geom_hline(yintercept=0,color="darkgrey") +
  geom_vline(xintercept=0,color="darkgrey") +
  labs(x = "Time", y = "Densities", color="Legend",title="Moose cap")+
  scale_color_manual(values = colors)

#wolf hunting 

sP <- 0.5
Kp <- 240
lP <- 2.68275
hp <- 300
sL <- 0.07
Kl <- 870
lL <- 2.68275
hl <- 400
sH <- 0.5
Kh <- 970
lH <- 3.1682
hH <- 1000
fC <- 1
eC <- 1.85
mC <- 0
fM <- 1.5
eM <- 0.6
b <- 0.8
g <- 0.38
dC <- 460
dM <- 46
mW <- 0.1


nsteps <- 200
pop.df.4 <- data.frame(time=1:nsteps,P=numeric(nsteps),L=numeric(nsteps),H=numeric(nsteps),C=numeric(nsteps),M=numeric(nsteps),W=numeric(nsteps))

pop.df.4 <- within(pop.df.4,{
  P[1] <- 240
  L[1] <- 870
  H[1] <- 970
  C[1] <- 7
  M[1] <- 25
  W[1] <- 8
})


nsteps <- 200
pop.df.4 <- data.frame(time=1:nsteps,P=numeric(nsteps),L=numeric(nsteps),H=numeric(nsteps),C=numeric(nsteps),M=numeric(nsteps),W=numeric(nsteps))

pop.df.4 <- within(pop.df.4,{
  P[1] <- 240
  L[1] <- 870
  H[1] <- 970
  C[1] <- 7
  M[1] <- 25
  W[1] <- 8
})


for(t in 2:nsteps){
  pop.df.4 <- within(pop.df.4,{
    ##Primary producers
    P.birth <- sP*P[t-1]*(1-P[t-1]/Kp)
    P.death <- lP*P[t-1]*C[t-1]/(hp+P[t-1])
    P[t] <- max(0,P[t-1] + P.birth - P.death) ##plants consumed by Caribou - the max is used to avoid negative values
    
    L.birth <- sL*L[t-1]*(1-L[t-1]/Kl)
    L.death <- lL*L[t-1]*C[t-1]/(hl+L[t-1])
    L[t] <- max(0,L[t-1] + L.birth - L.death) ##lichen consumed by Caribou
    
    P.birth <- sH*H[t-1]*(1-H[t-1]/Kh)
    P.death <- lH*H[t-1]*M[t-1]/(hH+H[t-1])
    H[t] <- max(0,H[t-1] + P.birth - P.death) ##plants consumed by Moose
    
    ##First trophic level
    C.growth <- C[t-1]*fC*(P[t-1]/(P[t-1]+hp))*(L[t-1]/(L[t-1]+hl))
    C.death <- C[t-1]*(1-P[t-1]/(P[t-1]+hp))*(1-L[t-1]/(L[t-1]+hl))
    C.pred <-  W[t-1]*eC*C[t-1]/(C[t-1]+dC) 
    C[t] <- max(0,C[t-1] + C.growth - C.death - C.pred)
    
    M.growth <- M[t-1]*fM*H[t-1]/(H[t-1]+hH)
    M.death <- M[t-1]*(1-H[t-1]/(H[t-1]+hH))
    M.pred <-  W[t-1]*eM*M[t-1]/(M[t-1]+dM)
    M[t] <- max(0,M[t-1] + M.growth - M.death - M.pred)
    
    
    ##Predator
    W.growth <- (b*(C[t-1]/(dC+C[t-1])) + b*(M[t-1]/(dM+M[t-1])))*W[t-1]
    W.death <- g*W[t-1]
    W[t] <- max(0,W[t-1] + W.growth - W.death - W[t-1]*mW)
  })
}


colors <- c("Plants"="brown","Lichen"="orange","Shrubs"="green","Caribou"="purple","Moose"="blue","Wolf"="red")
ggplot(data = pop.df.4)+
  geom_line(mapping=aes(x=time,y=P,color="Plants")) +
  geom_line(mapping=aes(x=time,y=L,color="Lichen")) +
  geom_line(mapping=aes(x=time,y=H,color="Shrubs")) +
  geom_line(mapping=aes(x=time,y=C,color="Caribou")) +
  geom_line(mapping=aes(x=time,y=M,color="Moose")) +
  geom_line(mapping=aes(x=time,y=W,color="Wolf")) +
  geom_hline(yintercept=0,color="darkgrey") +
  geom_vline(xintercept=0,color="darkgrey") +
  labs(x = "Time", y = "Densities", color="Legend",title="Wolf hunting")+
  scale_color_manual(values = colors)

#moose hunting 

sP <- 0.5
Kp <- 240
lP <- 2.68275
hp <- 300
sL <- 0.07
Kl <- 870
lL <- 2.68275
hl <- 400
sH <- 0.5
Kh <- 970
lH <- 3.1682
hH <- 1000
fC <- 1
eC <- 1.85
mC <- 0
fM <- 1.5
eM <- 0.6
b <- 0.8
g <- 0.38
dC <- 460
dM <- 46

##add hunting parameter
mM <- 0.1


nsteps <- 200
pop.df.5 <- data.frame(time=1:nsteps,P=numeric(nsteps),L=numeric(nsteps),H=numeric(nsteps),C=numeric(nsteps),M=numeric(nsteps),W=numeric(nsteps))

pop.df.5 <- within(pop.df.5,{
  P[1] <- 240
  L[1] <- 870
  H[1] <- 970
  C[1] <- 7
  M[1] <- 25
  W[1] <- 8
})


for(t in 2:nsteps){
  pop.df.5 <- within(pop.df.5,{
    ##Primary producers
    P.birth <- sP*P[t-1]*(1-P[t-1]/Kp)
    P.death <- lP*P[t-1]*C[t-1]/(hp+P[t-1])
    P[t] <- max(0,P[t-1] + P.birth - P.death) ##plants consumed by Caribou - the max is used to avoid negative values
    
    L.birth <- sL*L[t-1]*(1-L[t-1]/Kl)
    L.death <- lL*L[t-1]*C[t-1]/(hl+L[t-1])
    L[t] <- max(0,L[t-1] + L.birth - L.death) ##lichen consumed by Caribou
    
    P.birth <- sH*H[t-1]*(1-H[t-1]/Kh)
    P.death <- lH*H[t-1]*M[t-1]/(hH+H[t-1])
    H[t] <- max(0,H[t-1] + P.birth - P.death) ##plants consumed by Moose
    
    ##First trophic level
    C.growth <- C[t-1]*fC*(P[t-1]/(P[t-1]+hp))*(L[t-1]/(L[t-1]+hl))
    C.death <- C[t-1]*(1-P[t-1]/(P[t-1]+hp))*(1-L[t-1]/(L[t-1]+hl))
    C.pred <-  W[t-1]*eC*C[t-1]/(C[t-1]+dC) 
    C[t] <- max(0,C[t-1] + C.growth - C.death - C.pred)
    
    M.growth <- M[t-1]*fM*H[t-1]/(H[t-1]+hH)
    M.death <- M[t-1]*(1-H[t-1]/(H[t-1]+hH))
    M.pred <-  W[t-1]*eM*M[t-1]/(M[t-1]+dM)
    M[t] <- max(0,M[t-1] + M.growth - M.death - M.pred - M[t-1]*mM)
    
    ##Predator
    W.growth <- (b*(C[t-1]/(dC+C[t-1])) + b*(M[t-1]/(dM+M[t-1])))*W[t-1]
    W.death <- g*W[t-1]
    W[t] <- max(0,W[t-1] + W.growth - W.death)
  })
}


colors <- c("Plants"="brown","Lichen"="orange","Shrubs"="green","Caribou"="purple","Moose"="blue","Wolf"="red")
ggplot(data = pop.df.5)+
  geom_line(mapping=aes(x=time,y=P,color="Plants")) +
  geom_line(mapping=aes(x=time,y=L,color="Lichen")) +
  geom_line(mapping=aes(x=time,y=H,color="Shrubs")) +
  geom_line(mapping=aes(x=time,y=C,color="Caribou")) +
  geom_line(mapping=aes(x=time,y=M,color="Moose")) +
  geom_line(mapping=aes(x=time,y=W,color="Wolf")) +
  geom_hline(yintercept=0,color="darkgrey") +
  geom_vline(xintercept=0,color="darkgrey") +
  labs(x = "Time", y = "Densities", color="Legend",title="Moose hunting")+
  scale_color_manual(values = colors)















