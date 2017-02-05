
## Testing Bayesian Estimatation of Ccrit from hypothetical data of Ci and Fv'Fm ###
# R 3.2.1
# Updated: 02__2017. Jonathan R Plean.
# Reviewed by: No one

### Script depends on two models for finding Ccrit 
# both models based off Relialbe Estimation of biochemical paramters of C3 leaf...
              ###          Gu et al (2010) PCE Figs. 8,9,10

# model estimates Ccrit (transition point btw light and carbon limited photosynthesis)
# alpha (y intercept on Fv'Fm as Ci approaches 0)
# 2 beta values  (beta[1] slope prior to Ccrit, beta[2] slope after Crit)

# Sources: collaborator 
#


require("coda")
require("rjags")


setwd("~/Desktop")
source("Ccrit_model.R")


#### Set up for rjags #########
parameters = c("alpha", "beta","Ccrit","tau")### pars to be monitored
adaptSteps = 1000             # Number of steps to "tune" the samplers.
burnInSteps = 5000            # Number of steps to "burn-in" the samplers.
nChains = 4                   # Number of chains to run.
DICsteps= 20000                # Number of steps of sample DIC
numSavedSteps= 5000       # Total number of steps in chains to save.
thinSteps=20                   # Number of steps to "thin" (1=keep every step).
nPerChain = ceiling( ( numSavedSteps * thinSteps ) / nChains ) # Steps per chain.
###################################
####  Hypothetical DATA   #####
## 12 individuals 4 genotpyes
## six intercepts
Al<-rep(seq(0.42,0.52,by=0.02), 2)
## six inital slopes
b1<-rep(seq(8e-3,9e-3, length.out=6),2)
##slope after Ccrit
b2<-0
# 6 Ccrit's
Ccrit<-c(20,22,24,30,32,34, 40,42,44, 50,52,54)
Cdat<-c(seq(0,40,5),seq(50,150,25))


dat <- matrix(NA, ncol =12, nrow = 14)
for(i in 1:12){
  dat[,i]<-ifelse(Cdat<Ccrit[i], Cdat*b1[i]+Al[i], Cdat*b2+Ccrit[i]*b1[i]+Al[i])
}
### plotting a subset of hyothetical data
plot(Cdat,dat[,1], ylim=c(0.4,0.9))
Col<-seq(1,6,1)
for(i in 1:5){
points(Cdat,dat[,i+1],col=Col[i+1])
}
for(i in 1:6){
abline(v=Ccrit[i], col=Col[i])
}
##########################################

##########################################

##########################################
###    non-hierarchical model set & sampling

##########################################
setwd("~/Desktop")
source("Ccrit_model.R")
### data organized for passing to jags model
datalist<-vector("list", 12)
for(j in 1:12){
  datalist[[j]]<-list(N=length(Cdat),
                      CiP=Cdat,  Fv..Fm.=dat[,j])
}
### create models & burn in
models<-vector("list", 12)
for(j in 1:12){
  models[[j]] <- jags.model(textConnection(Ccrit_model), 
                            data = datalist[[j]], n.chains=nChains , n.adapt=adaptSteps)
  update(models[[j]], burnInSteps)
}
### sampling all models, calc sigma and convert to martrix 
mcmcsamples<-vector("list", nChains); 
for(j in 1:12){
  mcmcsamples[[j]] <- coda.samples(models[[j]],
                                   variable.names=parameters,
                                   n.iter=nPerChain , thin=thinSteps )
}
mcmcChain<-vector("list", 12); sigma<-vector("list", 12)
for(i in 1:12){
  mcmcChain[[i]] <- as.matrix( mcmcsamples[[i]])
  sigma[[i]] =1  / sqrt( mcmcChain[[i]][, "tau" ] )
  mcmcChain[[i]] = as.data.frame(cbind( mcmcChain[[i]], sigma[[i]] ))
}
# calc medians
meds<-matrix(, nrow = 12, ncol = 6)
for(j in 1:12){
  meds[j,]<- apply(mcmcChain[[j]], 2, median)
  rownames(meds)<-seq(1,12,1)
  colnames(meds)<-c("Ccrit","alpha","beta1","beta2","tau","sigma")  
}
meds<-as.data.frame(meds)
### quick comparison median posterior est against hypothetical values
plot(meds$Ccrit,Ccrit)
abline(0,1)

### convergence diagnostics
GD<-matrix(, nrow = 6, ncol=5)
for(j in 1:6){
  GD[j,]<- gelman.diag(mcmcsamples[[j]])$psrf[1]
  rownames(GD)<-seq(1,6,1)
  colnames(GD)<-c("Ccrit","alpha","beta1","beta2","tau")  
}
GDm<-matrix(, nrow = 6, ncol=1)
for(j in 1:6){
  GDm[j,]<- gelman.diag(mcmcsamples[[j]])$mpsrf[1]
  rownames(GDm)<-seq(1,6,1)
  colnames(GDm)<-c("multi_var")  
}
GD<-cbind(GD,GDm)
GD





###   hierarchical model set & sampling
setwd("~/Desktop")
source("Ccrit_model_hier.R")
datalisthier<-list(N=168, Ngeno=4,  CiP=rep(Cdat,12) ,
                   geno=c(rep(c(1:4), each=42)),
                     Fv..Fm.=c(dat[,1],dat[,2],dat[,3],dat[,4],dat[,5],dat[,6],
                               dat[,7],dat[,8],dat[,9],dat[,10],dat[,11],dat[,12]))
### create models, burn in & sample
model <- jags.model(textConnection(Ccrit_model_hier),
                    data = datalisthier, n.chains=nChains , n.adapt=adaptSteps)
update(model, burnInSteps)
mcmc_samples<- coda.samples(model,
                            variable.names=parameters,
                            n.iter=nPerChain , thin=thinSteps )

mcmcChain= as.matrix( mcmc_samples)
gelman.diag(mcmc_samples)
medsH<-as.data.frame(apply(mcmcChain,2,median))


### quick comparison median posterior est against hypothetical values
plot(Ccrit, rep(c(medsH[1:4,1]), each=3))
abline(0,1)






