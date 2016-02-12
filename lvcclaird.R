#this file is the mean field model in Laird et al 2009.
#JTB 256 (2009) 90-95

library(deSolve)
rm(list=ls())

#the size of the system is 
size<-5

#define alphas zeros indicate no interaction
alphas<-matrix(0,nrow=size,ncol=size) #competition coefs.

#define the state variables (populations) and their initial conditions are
#since this is a matrix notation system, use a named vector for 
#graphing purposes.
#state<-rep(0.2, size)
state<-c(0.1,0.2,0.2,0.2,0.3)
names(state)<-paste("N",1:size,sep="")

#to assign the competition coefficients, it was assumed 
#pairwise competitions in competitive exclussion
#then if A->B kB/kA < alpha21 and  kB/kA > alpha12
#        A<-B kB/kA > alpha21  and kA/kB < alpha12
#        ####       #                    #
##alpha[1,1]<- 0
alphas[1,2]<- 1 # 1->2
alphas[2,1]<- -1 ##1->2  
alphas[1,3]<- 1 # 1->3
alphas[3,1]<- -1 ##1->3
alphas[1,4]<- -1 # 1<-4
alphas[4,1]<- 1 ##1<-4
alphas[1,5]<- -1 # 1<-5
alphas[5,1]<- 1 ##1<-5
##  alpha[2,2]<- 0
alphas[2,3]<- 1 # 2->3
alphas[3,2]<- -1 ##2->3
alphas[2,4]<- 1 # 2->4
alphas[4,2]<- -1 ##2->4
alphas[2,5]<- -1 # 2<-5
alphas[5,2]<- 1 ##2<-5
##  a33<-0
alphas[3,4]<- 1 # 3->4
alphas[4,3]<- -1 ##3->4
alphas[3,5]<- 1 #3 ->5
alphas[5,3]<- -1 ##3->5
##  a44<-0
alphas[4,5]<- 1 # 4->5
alphas[5,4]<- -1 ##4<-5
#alphas[5,5]<-0


#gather all parameters in a list
parameters<-list(alphas=alphas)

#define the vector of output times
times<-seq(from=0,to=100,by=0.5)

#define the differential equation

lvccmat<-function(t,Y,parameters)
{
  y<- Y #transform Y in a local vector
  #use the local vector to construct the system of ODEs
  #and assign it to the vector of ODEs
  dY<- colSums(outer(y,y,FUN ="*")*alphas)
  #return the whole vector for the output
  return(list(dY))
}

#now do the simulation

out<-ode(y=state,times=times,func=lvccmat)
#this plots the output
plot(out)
