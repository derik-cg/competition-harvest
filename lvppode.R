#This file simulate a 2 population competition with lotka-Volterra
#predator-prey equations. This is a base for constructing an n-population
#model

library(deSolve)

# define the parameter vector
parameters<-c(r1=1.2,
              r2=1.4,
              a12=0.6,
              a21=0.7)

#define the state variables (populations) and their initial conditions are
state<-c(V=5, P=5)

#define the vector of times for the output
times<-seq(from=0,to=100,by=0.5)

# define the equations
lvppode<-function(t,state,parameters)
{
  with(as.list(c(state,parameters)),
       {
         #this defines the differential equations
         dV<-  r1*V - (a12*V*P)
         dP<- -r2*P + (a12*V*P)
         #this defines the output
         list(c(dV,dP))
       }) #this parenthesis ends the with
}

out<-ode(y=state, times=times, func=lvppode, parms=parameters)
