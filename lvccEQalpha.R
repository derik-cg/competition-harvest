#this files simulates a community of five populations competing
#according to the Lotka-Volterra assumptions.
#the objective is to try intransitive topologies of the
#competition network
#this applies laird et al 2009 on a lv system

#the objective of this file is to examine the stability of the
#equilibrium point when alphas are changed one at a time.

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
#do a for loop for one of the alpas
#define first a vector of alphas to simulate from
vecalphas<-seq(from=1.1,to=10.9,length.out = 20)
vbig<-NULL
vsma<-NULL
vtf<-NULL
for(sma in vecalphas)
{
  big<-1.5
  #sma<-0.5
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
  
  equilibria<-function(df)
  {
    sz<-dim(df)
    #this function takes the output of an ode solve, clips out
    #the last 50 entries and find if they are at equilibria
    eq<-out[(sz[1]-50):sz[1],2:sz[2]]
    #check if the equilibrium was reached  
    minmax<-function(x){c(min=min(x),max=max(x))}
    eqminmax<-apply(eq,2,minmax)
    abs(eqminmax["min",]-eqminmax["max",])<0.01
  }
  vbig<-c(vbig,big)
  vsma<-c(vsma,sma)
  vtf<-rbind(vtf,equilibria(out))
}
dfequi<-data.frame(alphaBIG=vbig,alphaSMALL=vsma,vequi=vtf)
