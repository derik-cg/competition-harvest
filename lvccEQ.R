#this files simulates a community of five populations competing
#according to the Lotka-Volterra assumptions.
#the objective is to try intransitive topologies of the
#competition network
#this applies laird et al 2009 on a lv system

#The objective of this file is to show that symmetric deviations
#in alphas from the ratio k1/k2 and k2/k1 preserve the numerical
#value of the equilibrium.
 
library(deSolve)
rm(list=ls())

#the size of the system is 
size<-5

# define the parameter vectors
rs<-rep(0.5,size) #growth rates
ks<-rep(100,size) #carrying capacities

#define alphas zeros indicate no interaction
alphas<-matrix(0,nrow=size,ncol=size) #competition coefs.

#define the state variables (populations) and their initial conditions are
#since this is a matrix notation system, use a named vector for 
#graphing purposes.
state<-rep(c(20,20), length.out=size)
names(state)<-paste("N",1:size,sep="")

#to assign the competition coefficients, it was assumed 
#pairwise competitions in competitive exclussion
#then if A->B kB/kA < alpha21 and  kB/kA > alpha12
#        A<-B kB/kA > alpha21  and kA/kB < alpha12
#        ####       #                    #
##alpha[1,1]<- 0
dist<-0.9
big<-1+dist
sma<-1-dist
alphas[1,2]<- big # 1->2
alphas[2,1]<- sma ##1->2  
alphas[1,3]<- big # 1->3
alphas[3,1]<- sma ##1->3
alphas[1,4]<- sma # 1<-4
alphas[4,1]<- big ##1<-4
alphas[1,5]<- sma # 1<-5
alphas[5,1]<- big ##1<-5
##  alpha[2,2]<- 0
alphas[2,3]<- big # 2->3
alphas[3,2]<- sma ##2->3
alphas[2,4]<- big # 2->4
alphas[4,2]<- sma ##2->4
alphas[2,5]<- sma # 2<-5
alphas[5,2]<- big ##2<-5
##  a33<-0
alphas[3,4]<- big # 3->4
alphas[4,3]<- sma ##3->4
alphas[3,5]<- big #3 ->5
alphas[5,3]<- sma ##3->5
##  a44<-0
alphas[4,5]<- big # 4->5
alphas[5,4]<- sma ##4<-5
#alphas[5,5]<-0


#gather all parameters in a list
parameters<-list(rs=rs,
                 ks=ks,
                 alphas=alphas)

#define the vector of output times
times<-seq(from=0,to=1000,by=0.5)

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
