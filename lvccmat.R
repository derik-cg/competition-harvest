#this files simulates an arbrtrary (>2) number of competing populations
# using matrix notation for the differential equations. No harvest

library(deSolve)
rm(list=ls())
# define the parameter vector

rs<-c(0.2,0.5) #growth rates
ks<-c(100,100) #carrying capacities
alphas<-matrix(0,nrow=2,ncol=2,byrow = T) #competition coefs.
alphas[1,2]<-0.5
alphas[2,1]<-1.5
#the size of the system is 
size<-length(rs)

#define the state variables (populations) and their initial conditions are
#since this is a matrix notation system, no not use names
state<-rep(10, size)

#gather all parameters in a list
parameters<-list(rs=rs,
              ks=ks,
              alphas=alphas,
              size=size)

#define the vector of output times
times<-seq(from=0,to=100,by=0.5)

#define the differential equation

lvccmat<-function(t,Y,parameters)
{
  y<- Y #transform Y in a local vector
  #use the local vector to construct the system of ODEs
  #and assign it to the vector of ODEs
  dY<- rs*y - (rs*y^2 + rs*colSums(outer(y,y,FUN ="*")*alphas))*(1/ks)
  #return the whole vector for the output
  return(list(dY))
}

#now do the simulation

out<-ode(y=state,times=times,func=lvccmat)
#this plots the output
plot(out)
