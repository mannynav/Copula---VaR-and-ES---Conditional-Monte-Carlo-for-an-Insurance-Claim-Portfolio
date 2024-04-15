
#Packages.
library(copula)
library(MASS)
library(QRM)
library(fitdistrplus)
library(CASdatasets)
library(pscl)
library(VineCopula)

#Get the data.
data(ukaggclaim)

#Claim severity data
ukclaim <- ukaggclaim$AvCost
summary(ukclaim)
skewness(ukclaim)
kurtosis(ukclaim)
ukclaim <- as.numeric(na.omit(ukclaim))

#Claim frequency data
ukfreq <- ukaggclaim$NClaims
var(ukfreq)
summary(ukfreq)
skewness(ukfreq)
kurtosis(ukfreq)
ukfreq <- ukfreq[ukfreq>0]

#Histogram of data
par(mfrow=c(1,2))
hist(ukfreq, main = c("Histogram of Claim Frequency"),xlab = c("Claims"))
hist(ukclaim,main = c("Histogram of Severities"),xlab = c("Severity (pounds)"))
par(mfrow=c(1,1))

#Fit for claim frequency.
fit_pois <- fitdist(ukfreq,"pois")
fit_nbinom <- fitdist(ukfreq, distr = "nbinom")
fit_geometric <- fitdist(ukfreq, distr = "geom")

#Summaries for claim frequency.
summary(fit_pois)
summary(fit_nbinom)
summary(fit_geometric)

#Goodness-of-fit tests for claim frequency.
gofstat(list(fit_pois,fit_nbinom,fit_geometric))

#Generalized linear model.
glmNB <- glm.nb(NClaims ~1, data = ukaggclaim)
glmNB

#We get the same mu as the fit dist function.
exp(glmNB$coefficients)

#Zero inflated model with negative binomial distribution
nbZI <- zeroinfl(NClaims ~ 1,dist = "negbin",data = ukaggclaim)
summary(nbZI)

#Mu
exp(4.2465)

#Size
exp(-0.5121)

#Compare fits using GLM and zero inflated model.
AIC(glmNB,nbZI)

#Fit for claim severity. This is be the marginal distribution for X_1,X_2,...,X_n.
fit_pareto <- fitdist(ukclaim, distr = "pareto", start = list(shape = 5, scale = 1))
fit_burr <- fitdist(ukclaim,distr = "burr",start = list(shape1 = 1, shape2 = 1, rate = 1))
fit_lognorm <- fitdist(ukclaim,distr = "lnorm")
fitdistr_weibull <- fitdist(ukclaim,"weibull",lower = c(0,0), start = list(shape = 0.05, scale = 0.1))
fit_gamma <- fitdist(ukclaim, distr = "gamma")
fit_exp <- fitdist(ukclaim, distr = "exp")

#Summaries for claim severity.
summary(fit_pareto)
summary(fit_burr)
summary(fit_lognorm)
summary(fitdistr_weibull)
summary(fit_gamma)
summary(fit_exp)

#Goodness-of-fit tests for severity.
gofstat(list(fit_pareto,fit_lognorm,fit_burr,fitdistr_weibull,fit_gamma,fit_exp))

###### Copula modelling #####

#Get empirical distribution for claim frequency and claim severity.
claimCDF<- ecdf(ukclaim)(ukclaim)
freqCDF<-ecdf(ukfreq)(ukfreq)

#Scatter plot for ECDF.
ret <- cbind(freqCDF,claimCDF)
plot(ret, xlab = "ECDF of Claim Frequency", ylab = "ECDF of Claim Severity", main = "Claim Frequency and Claim Severity ECDF Scatterplot")

#Fit the ECDF of the returns to bivariate copulas using BiCopSelect function.
?BiCopSelect

#Switch families to get appropriate copula.
CopulaFit<-BiCopSelect(freqCDF,claimCDF, family = 7,rotations = FALSE,indeptest = TRUE)
summary(CopulaFit)

