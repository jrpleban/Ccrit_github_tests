###############################
###      Bayesian Model    ####
### For Ccrit model ####
####
######  & no temp dependency ###
###############################
Ccrit_model <-
"model {
    for (i in 1:N){
        
           Fv..Fm.[i] ~ dnorm( mu.Fv..Fm.[i] , tau )
           mu.Fv..Fm.[i] <- alpha + beta[J[i]] * (CiP[i] - Ccrit)
           J[i] <- 1 + step(CiP[i] - Ccrit)
    }
        tau ~ dgamma(0.001, 0.001)
        alpha ~ dnorm(0.0, 1.0E-6)
        beta[1]~ dnorm(0.0,1.0E-6) #Uninformative
        beta[2]~ dnorm(0.0,1.0E-6) #Uninformative
        Ccrit ~ dnorm(30,.01)
}

" 