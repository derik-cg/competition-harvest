#this files simulates a community of five populations competing
#according to the Lotka-Volterra assumptions.
#the objective is to try intransitive topologies of the
#competition network
#this applies laird et al 2009 on a lv system

#The objective of this file is to explore ASYMETTRIC deviations
#in alphas from the ratio k1/k2 and k2/k1 in relation to the numerical
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

#define the differential equaiton
lvccmat<-function(t,Y,parameters)
{
  y<- Y #transform Y in a local vector
  #use the local vector to construct the system of ODEs
  #and assign it to the vector of ODEs
  dY<- rs*y - (rs*y^2 + rs*colSums(outer(y,y,FUN ="*")*alphas))*(1/ks)
  #return the whole vector for the output
  return(list(dY))
}
#define function that measures amplitude of oscillations
amplitude<-function(mat)
{
  mx<-apply(mat,2,max)
  mn<-apply(mat,2,min)
  dife<-mx-mn
  return(dife)
}

#define the vector of output times
times<-seq(from=0,to=1000,by=0.5)

#initialize output of oscilations
out.osc<-c(dist=NA,eq1=NA,eq2=NA,eq3=NA,eq4=NA,eq5=NA,N1=NA,N2=NA,N3=NA,N4=NA,N5=NA)
## start a loop for the symmetrical differences
distseq<-seq(from=0.0001, to=0.9999, length.out = 50)
for (i in distseq)
{
  #to assign the competition coefficients, it was assumed 
  #pairwise competitions in competitive exclussion
  #then if A->B kB/kA < alpha21 and  kB/kA > alpha12
  #        A<-B kB/kA > alpha21  and kA/kB < alpha12
  #        ####       #                    #
  ##alpha[1,1]<- 0
  big<-1+i
  sma<-1-0.5
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
  
  
  #now do the simulation
  
  out<-ode(y=state,times=times,func=lvccmat)
  #this plots the output
  #plot(out)
  # take the last 50 rows
  last<-out[1550:2000,2:6]
  #measure for oscillations
  ampl<-amplitude(last)
  #estimate equilibria with average
  avg<-apply(last,2,mean)
  #name the entries of the vector#
  names(avg)<-c("eq1","eq2","eq3","eq4","eq5")
  #associate dist, equilibria and amplitude
  osc<-c(dist=i,avg,ampl)
  
  out.osc<-rbind(out.osc,osc)
}
#eliminate fist row
num<-dim(out.osc)
out.osc<-out.osc[2:(num[1]),]
#this plot relates the distance to the equilibrium point.
plot((out.osc[,1]),out.osc[,2])
