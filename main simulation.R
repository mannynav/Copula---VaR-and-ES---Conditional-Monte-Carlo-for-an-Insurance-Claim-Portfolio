
#Packages.
library(actuar)
library(pracma)

#Number of simulations 250,000.
R<- 250000

#Create vector to store sum of X1+X2+...+X_{n-1}.
vectS <- c(0)

for (i in 1:R) {
  #Claim frequency distribution.
  N <- rnbinom(1,size = fit_nbinom$estimate[1], mu = fit_nbinom$estimate[2])
  
  #Generate V1,V2,...Vn
  V <- rexp(N,1)
  
  #Set seed for replication.
  set.seed(i)
  
  #Random variable Z >0.
  Z <- rexp(1,1)
  #Zclay <- rgamma(1,1/0.55)
  
  #V divided by Z.
  VdZ <- V/Z
  
  #Generator functions for Archimedean copulas.
  #psiClay <- 1/((1+VdZ*0.55)^(1/0.55))
  #psiBB7 <- 1 - (1 - (1/(VdZ+1))^(1/0.55)) ^ (1/1)
  psiBB1 <- (1/(VdZ^(1/1)+1))^(1/0.55)
  #psiJoe <- 1 - (1-exp(-VdZ))^(1/1)
  
  #Get inverse of marginal distribution F.
  Finvs <- 1 - as.vector(pburr(psiBB1,fit_burr$estimate[1],fit_burr$estimate[2],fit_burr$estimate[3]))
  
  #Vector to store n-1 values of Finvs.
  vectS[i] <- sum(Finvs[1:length(Finvs)-1])
}

#Parameters for BB1 copula.
thetaBB1<-.55
betaBB1 <- 1

#Parameters for BB7 copula.
thetaBB7<-1
betaBB7 <- 0.55

#Parameter for Clayton copula.
thetaClay<-0.55

#Parameter for Joe copula.
thetaJoe <-1

#Estimated CDF.
Fncond <- function(x,vectS) {
  s<-c(0) 
  for (i in 1:R) {
    set.seed(i)
    Z <- rexp(1,1)
    
    #F(x-S_{n-1})
    inputF <- pburr(x-vectS[i],fit_burr$estimate[1],fit_burr$estimate[2],fit_burr$estimate[3])
    
    #Inverse generators for each copula.
    invGeneratorBB1 <- (inputF^(-thetaBB1)-1)^(betaBB1)
    invGeneratorBB7 <- ((1-(1-inputF)^thetaBB7)^(-betaBB7)) - 1
    invGeneratorClay <- (1/((inputF)^(thetaClay)) - 1)/thetaClay
    invGeneratorJoe <- (-log(1-(1-inputF)^(thetaJoe)))
    s[i] <- exp(-Z*invGeneratorBB1)
  }
  return(mean(s))
}

#Compute the inverse of Fncond.
Fncondinvs <- function(x) {
  sINV<-c(0) 
  for (i in 1:R) {
    set.seed(i)
    Z <- rexp(1,1)
    
    #F(x-S_{n-1})
    inputF <- pburr(x-vectS[i],fit_burr$estimate[1],fit_burr$estimate[2],fit_burr$estimate[3])
    
    #Inverse generators for each copula.
    invGeneratorBB1 <- (inputF^(-thetaBB1)-1)^(betaBB1)
    invGeneratorBB7 <- ((1-(1-inputF)^thetaBB7)^(-betaBB7)) - 1
    invGeneratorClay <- (1/((inputF)^(thetaClay)) - 1)/thetaClay
    invGeneratorJoe <- (-log(1-(1-inputF)^(thetaJoe)))
    sINV[i] <- exp(-Z*invGeneratorBB1)
  }
  return(1-mean(sINV))
}

#Vectorize inverse Fncondinvs
FninvscondVect <- Vectorize(Fncondinvs, "x")

#Quantile.
alpha <- 0.95

#Solve Fncond(1/(1-q)-1,vectS,theta) = 0.95 for q, store as qu
qu <- uniroot(f = function(q) alpha - Fncond(1/(1 - q) - 1,vectS), interval = c(0,1),check.conv = TRUE)$root

#This is the VaR estimate
VaR <- 1 /(1 - qu) - 1

#Integrate the inverse distribution function from the VaR estimate to infinity to get the estimate for m.
fun <- function(x) FninvscondVect(x)
m <- integral(fun,VaR,Inf,reltol = 1e-8)

#Compute expected shortfall.
ES <- VaR + (m)/(1-alpha)

