#this file simulates the lotka volterra competition equations without 
#harvest
library(deSolve)
rm(list=ls())
# define the parameter vector
parameters<-c(r1=1.2,
              r2=1.4,
              k1=100,
              k2=90,
              a12=0.6,
              a21=0.7)

#define the state variables (populations) and their initial conditions are
state<-c(N1=5, N2=5)

#define the vector of times for the output
times<-seq(from=0,to=100,by=0.5)

# define the equations
lvcc<-function(t,state,parameters)
{
  with(as.list(c(state,parameters)),
       {
         #this defines the differential equations
         dN1<- r1*N1*(k1 - N1 - a12*N2)/k1
         dN2<- r2*N2*(k2 - N2 - a21*N1)/k2
         #this defines the output
         list(c(dN1,dN2))
       }) #this parenthesis ends the with
}

out<-ode(y=state, times=times, func=lvcc, parms=parameters)
