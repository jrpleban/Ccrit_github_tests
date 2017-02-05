###############################
###  Hierarchical Bayesian Model    ####
###        For Ccrit model          ####
###############################
Ccrit_model_hier <-
  "model {
for (i in 1:N){

Fv..Fm.[i] ~ dnorm( mu.Fv..Fm.[i] , tau )
mu.Fv..Fm.[i] <- alpha[geno[i]] + beta[J[i],geno[i]] * (CiP[i] - Ccrit[geno[i]])
J[i] <- 1 + step(CiP[i] - Ccrit[geno[i]])
}
for (g in 1:Ngeno){
beta[1,g]~ dnorm(mu.beta, tau.beta)
beta[2,g]~ dnorm(mu.beta, tau.beta)
alpha[g] ~ dnorm(mu.alpha, tau.alpha)T(0,)
Ccrit[g] ~ dnorm(mu.Ccrit,tau.Ccrit)T(0,)
}

mu.alpha ~ dnorm(1, 0.1)T(0,)
mu.beta~ dnorm(0,0.01) #Uninformative
mu.Ccrit ~ dnorm(30, 0.01)T(0,)
tau.alpha ~ dnorm(.2, 1)
tau.beta~ dnorm(1.0,1) #Uninformative
tau.Ccrit ~ dnorm(4,.25)
tau ~ dgamma(0.001, 0.001)
}


" 