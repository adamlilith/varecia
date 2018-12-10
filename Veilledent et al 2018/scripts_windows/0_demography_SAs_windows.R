#===========================================================
# R script: Demographic model for the phcf study areas
#
# Copyright (C) March 2012:
# Ghislain Vieilledent <ghislain.vieilledent@cirad.fr>
# Clovis Grinand <clovis@goodplanet.org>
# Oliver Gardi <oliver.gardi@helvetas.org>
# License: GPL 3 (see license_windows.txt)
#
# The following code is distributed in the hope that it will be
# useful, but WITHOUT ANY WARRANTY.
#
# Reference: Vieilledent G., C. Grinand and
# R. Vaudry. 2013. Forecasting deforestation and carbon emissions in
# tropical developing countries facing demographic expansion: a case
# study in Madagascar. Ecology and Evolution. DOI: 10.1002/ece3.550
#
# This R script describes the steps to follow in order to estimate the
# parameters of the demographic model using the demography() function
# of the phcfM R package.
#
#==========================================================

#= Set working directory
# setwd("/path/to/0_data_demogaphy_SAs_windows.R")

#= Libraries
library(phcfM) # For demography() function

#= Import data
dataf <- read.table(file="./data/0_data_demography_SAs.txt",sep="\t",header=TRUE)

#= Create a time variable with t=0 for date 1993
dataf$time <- dataf$date-1993

#= Firaisana category
Levels.Fir <- sort(unique(dataf$codefir))
nfir <- length(unique(dataf$codefir)) # 253

#=================
#= Spaghetti plots
pdf(file="./outputs/demography/population_evolution.pdf")
par(cex=1.2,mar=c(5,5,2,1))
plot(x=dataf$time,
     y=dataf$pop,
     xlim=c(0,16),
     ylim=c(0,100000),
     type=c("p"),
     col="white",
     xlab="Year",
     ylab="Population of the Firaisana\n(in thousand)",
     axes=FALSE
     )
axis(side=1,at=c(0,11:16),labels=c(1993,2004,2005,2006,2007,2008,2009),las=3,tcl=0.5)
axis(side=2,at=seq(from=0,to=100000,by=20000),labels=as.character(seq(from=0,to=100,by=20)),tcl=0.5)
#= Lines
for (i in 1:nfir) {
  lines(x=dataf$time[dataf$codefir==Levels.Fir[i]],
        y=dataf$pop[dataf$codefir==Levels.Fir[i]],
        col=grey(0.5),
        lty=1,
        lwd=1)
}
#= Points
points(dataf$time,dataf$pop,col="black",pch=16)
dev.off()

#= Exponential demographic model
# dP/dt = theta*P, with theta equal to the "mean population growth rate"
# P(t) = rho*exp(theta*t)
# log(P_it) = (beta_0+b0_i) + (theta+b1_i)*t + eps_it
dataf$lP <- log(dataf$pop)
dataf$fir <- as.factor(as.character(dataf$codefir))

mod <- demography(fixed=lP~time,random=~1,group="fir",data=dataf,seed=1234,r=1,R=diag(1,1),burnin=1000,mcmc=1000,thin=1)
## plot(dataf$lP,mod$Y.pred)
rho.hat <- mean(mod$mcmc[,colnames(mod$mcmc)=="beta.time"])  # population growth rate: ~3.39%.yr-1

#= Output for population growth rate
sink("./outputs/demography/pop_growth_rate.txt")
cat("population growth rate (%.yr-1):",rho.hat)
sink()
#= Graphics for MCMC and posterior distributions
pdf("./outputs/demography/MCMC_posterior.pdf")
plot(mod$mcmc)
dev.off()

#==========================
#= Predictions by Firaisana
pdf(file="./outputs/demography/predictions.pdf")
par(cex=1.2,mar=c(5,5,2,1))
plot(x=dataf$time,
     y=dataf$pop,
     xlim=c(0,16),
     ylim=c(0,100000),
     type=c("p"),
     col="white",
     xlab="Year",
     ylab="Population of the Firaisana\n(in thousand)",
     axes=FALSE
     )
axis(side=1,at=c(0,11:16),labels=c(1993,2004,2005,2006,2007,2008,2009),las=3,tcl=0.5)
axis(side=2,at=seq(from=0,to=100000,by=20000),labels=as.character(seq(from=0,to=100,by=20)),tcl=0.5)
#= Lines
for (i in 1:nfir) {
  lines(x=dataf$time[dataf$codefir==Levels.Fir[i]],
        y=dataf$pop[dataf$codefir==Levels.Fir[i]],
        col=grey(0.7),
        lty=1,
        lwd=1)
}
#= Points
points(dataf$time,dataf$pop,col=grey(0.7),pch=16)

#= Predictions for each Firaisana
for (i in 1:nfir) {
  # x
  t.seq <- seq(from=0,to=max(dataf$time[dataf$codefir==Levels.Fir[i]]),length.out=100)
  # sigma2
  sigma2.hat <- mean(mod$mcmc[,colnames(mod$mcmc)=="sigma2"])
  # alpha0
  beta0.hat <- mean(mod$mcmc[,colnames(mod$mcmc)=="beta.(Intercept)"])
  b0.L.Fir <- paste("b.(Intercept).",Levels.Fir[i],sep="")
  b0.hat <- mean(mod$mcmc[,colnames(mod$mcmc)==b0.L.Fir])
  alpha0.hat <- beta0.hat+b0.hat+(sigma2.hat)/2 # Back transforming log-log relationship implies sigma.hat/2 for alpha0.hat !!
  # theta
  beta1.hat <- mean(mod$mcmc[,colnames(mod$mcmc)=="beta.time"])
  alpha1.hat <- beta1.hat
  # Predictions
  y.seq <- exp(alpha0.hat+alpha1.hat*t.seq)
  lines(t.seq,y.seq,col="black",lty=1,lwd=1)
}
dev.off()

#====================================================================
# Mean predictions with confidence envelops for Firaisana variability
pdf(file="./outputs/demography/mean_predictions.pdf")
par(cex=1.2,mar=c(5,5,2,1))
plot(x=dataf$time,
     y=dataf$pop,
     xlim=c(0,16),
     ylim=c(0,40000),
     type=c("p"),
     col="white",
     xlab="Year",
     ylab="Population of the Firaisana\n(in thousand)",
     axes=FALSE
     )
axis(side=1,at=c(0,11:16),labels=c(1993,2004,2005,2006,2007,2008,2009),las=3,tcl=0.5)
axis(side=2,at=seq(from=0,to=40000,by=10000),labels=as.character(seq(from=0,to=40,by=10)),tcl=0.5)
# Lines
for (i in 1:nfir) {
  lines(x=dataf$time[dataf$codefir==Levels.Fir[i]],
        y=dataf$pop[dataf$codefir==Levels.Fir[i]],
        col=grey(0.7),
        lty=1,
        lwd=1)
}
# Points
points(dataf$time,dataf$pop,col=grey(0.7),pch=16)

#===========================
#= Parameter for predictions

#= beta: fixed effects
beta0.mcmc <- mod$mcmc[,colnames(mod$mcmc)=="beta.(Intercept)"]
beta1.mcmc <- mod$mcmc[,colnames(mod$mcmc)=="beta.time"]
#= residual errors
sigma2.mcmc <- mod$mcmc[,colnames(mod$mcmc)=="sigma2"]
epsilon.mcmc <- rnorm(1000,0,sqrt(sigma2.mcmc)) 
#= Variance for random effects
Vb.mcmc <- mod$mcmc[,colnames(mod$mcmc)=="VCV.(Intercept).(Intercept)"]
#= b0: random effects
b0.mcmc <- rnorm(1000,mean=0,sqrt(Vb.mcmc))

#======================
#= Data for predictions
t.seq <- seq(from=0,to=max(dataf$time),length.out=100)
y.pred <- matrix(0,nrow=1000,ncol=length(t.seq))
for (i in 1:length(t.seq)) {
  y.pred[,i] <- exp(beta0.mcmc+b0.mcmc+beta1.mcmc*t.seq[i]+epsilon.mcmc)
}
#= quantiles
y.seq.q025 <- apply(y.pred,2,quantile,0.025)
y.seq.q975 <- apply(y.pred,2,quantile,0.975)
y.seq.mean <- apply(y.pred,2,mean)
#= lines
lines(t.seq,y.seq.q025,col="black",lty=2,lwd=2)
lines(t.seq,y.seq.q975,col="black",lty=2,lwd=2)
lines(t.seq,y.seq.mean,col="black",lty=1,lwd=2)

dev.off()

#=========================================================
# Cleaning folder

# Uncomment the following line if you want to clean the output folder

# shell("rm -R ./outputs/demography/*") 

#= End of script
#==========================================================================================#

