#this file is the differential equations for the lotka-volterra 
#competition with two populations and harvest

library(deSolve)
rm(list=ls())
# define the parameter vector
parameters<-c(r1=1.2,
              r2=1.4,
              a12=0.6,
              a21=0.7,
              h1=0.01,
              h2=0.01)

#define the state variables (populations) and their initial conditions are
state<-c(P=5, V=5)

#define the vector of times for the output
times<-seq(from=0,to=100,by=0.5)

# define the equations
lvpphv<-function(t,state,parameters)
{
  with(as.list(c(state,parameters)),
       {
         #this defines the differential equations
         dV<-  r1*V - (a12*V*P) - h1
         dP<- -r2*P + (a12*V*P) - h2
         #this defines the output
         list(c(dV,dP))
       }) #this parenthesis ends the with
}

out<-ode(y=state, times=times, func=lvpphv, parms=parameters)
