#==========================================================
# R script: Deforestation model.
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
# This R script describes the steps to follow in order to fit the
# deforestation model using the deforestation() function of the phcfM
# R package.
#
#==========================================================

#= Set working directory
# setwd("/path/to/2_deforestation_model_SA4_windows.R")

#= Library
library(phcfM)
library(vcd) # for Kappa

## #= Directories
## system("rm -r ./outputs/model")
## system("mkdir ./outputs/model")

#= Importing data on deforestation
data <- read.table("./data/1_data_deforestation_SA4.txt", sep="\t", header=TRUE)
colnames(data)
summary(data)

#==============================================
# Computing the population at the starting date
# Model: P(t)=K*exp(rho*t)

#= K, population at t=1993
rho <- 0.03389765 # Population growth rate
data$K_pop <- data$fkt_pop/exp(rho*(data$fkt_date-1993))

#= Population and density at the starting dates
num2date <- function(x) {
  x <- as.character(x)
  a <- substring(x,1,4)
  b <- substring(x,5,6)
  c <- substring(x,7,8)
  paste(a,b,c,sep="-")
}
data$date_image <- num2date(data$date_image)
data$diff.time <- as.numeric(difftime(data$date_image,"1993-06-01",units="days")/365.25)
data$pop <- data$K_pop*exp(rho*(data$diff.time))
data$dens <- data$pop/data$fkt_area

#= Remove extreme values of dens which can act as bad points lever for the regression
q99.dens <- quantile(data$dens,0.99)
data <- data[data$dens<=q99.dens,]

#======================================================
#= Mean deforestation rates on the data-set

model.mean <- deforestation(formula=statut_tp1~1, data=data, interval=data$interval,
                            verbose=1, seed=1234, tune=1, burnin=1000, mcmc=1000, thin=1)
MCMC.mean <- as.matrix(model.mean$mcmc)
theta.mean <- mean(inv.logit(MCMC.mean[,1])) # posterior mean of the annual deforestation rate
theta.prim.mean <- sum(data$statut_tp1==1)/(sum(data$statut_tp1==1)+sum(data$statut_tp1==0)) # ML estimates of the deforestation rate
thetas <- data.frame(theta.mean=theta.mean,theta.prim.mean=theta.prim.mean,row.names=NULL)

write.table(thetas,"./outputs/model/theta_thetaprim.txt",sep="\t",row.names=FALSE,quote=FALSE)

#======================================================
#= Mean deforestation rates evolution on the two periods of time

#= Selecting data for the 1st period
data.1 <- data[data$period==1,]
model.1 <- deforestation(formula=statut_tp1~1, data=data.1, interval=data.1$interval,
                            verbose=1, seed=1234, tune=1, burnin=1000, mcmc=1000, thin=1)
MCMC.1 <- as.matrix(model.1$mcmc)
theta.1 <- mean(inv.logit(MCMC.1)) # posterior mean of the annual deforestation rate

#= Selecting data for the 2sd period
data.2 <- data[data$period==2,]
model.2 <- deforestation(formula=statut_tp1~1, data=data.2, interval=data.2$interval,
                            verbose=1, seed=1234, tune=1, burnin=1000, mcmc=1000, thin=1)
MCMC.2 <- as.matrix(model.2$mcmc)
theta.2 <- mean(inv.logit(MCMC.2[,1])) # posterior mean of the annual deforestation rate

theta.period <- data.frame(period1=theta.1,period2=theta.2,row.names=NULL)

write.table(theta.period,"./outputs/model/theta_period.txt",sep="\t",row.names=FALSE,quote=FALSE)

#======================================================
#= Effect of the population

model.intensity <- deforestation(formula=statut_tp1~dens, data=data, interval=data$interval,
                            verbose=1, seed=1234, tune=1, burnin=1000, mcmc=1000, thin=1)

summary(model.intensity$mcmc)

#= Parameter estimates
MCMC.intensity <- as.matrix(model.intensity$mcmc)
beta.hat <- apply(MCMC.intensity,2,mean)
sd.hat <- apply(MCMC.intensity,2,sd)
q25.hat <- apply(MCMC.intensity,2,quantile,0.025)
q75.hat <- apply(MCMC.intensity,2,quantile,0.975)
sign0 <- rep(0,length(q75.hat))
sign0[q25.hat*q75.hat>0] <- 1
Estimates.model.intensity <- data.frame(Parameter=colnames(MCMC.intensity),beta=beta.hat,sd=sd.hat,sign=sign0,row.names=NULL)

write.table(Estimates.model.intensity,"./outputs/model/beta.txt",sep="\t",row.names=FALSE,quote=FALSE)

#================================
#= Modelling approach
#
# Une réalisation de la déforestation à t2:       Y_i = B(Theta_i')
# Probabilité du Pixel i d'être déforesté à t2:   Theta_i' = 1 - (1 - Theta_i)^T_i
# logit Model taux déforestation anuelle:         logit(Theta_i) = B0 + B1X1 + ... BnXn
#

#= Selecting data for the second period only
data <- data[data$period==2,]

#= Number of cross-validation
ncross <- 10
f <- list()

#= Variable names
varnames <- c("as.factor(prot_area)","altitude",
              "dist_road","dist_town","dist_patch","dist_fedge",
              "as.factor(fragindex)","dens","fkt_poors","fkt_cattle","fkt_mines","fkt_fert")
nvar <- length(varnames)

#= Formulas
for (v in 1:nvar) {
  f[[v]] <- formula(paste("statut_tp1~",paste(varnames[-c(v)],collapse="+"),sep="")) # we suppress 1 variable
}
f[[nvar+1]] <- formula(paste("statut_tp1~",paste(varnames,collapse="+"),sep="")) # full model
f[[nvar+2]] <- formula("statut_tp1~1") # null model
nmodels <- nvar+2

#=============================
#= Model with 100% of the data

#= Indexes
Deviance <- array(NA,dim=c(nmodels))

#= Loops on models
for (m in 1:nmodels) {                       

  #= All the data
  model.spatial <- deforestation(formula=f[[m]], data=data, interval=data$interval,
                                 verbose=1, seed=1234, tune=1, burnin=1000, mcmc=1000, thin=1)
  MCMC <- as.matrix(model.spatial$mcmc)

  #= Parameter estimates
  gamma.hat <- apply(MCMC,2,mean)

  #= Predictions on test data
  mf.fixed <- model.frame(formula=f[[m]],data=data)
  X <- model.matrix(attr(mf.fixed,"terms"),data=mf.fixed)
  theta.hat <- inv.logit(X%*%gamma.hat)
  theta.prim <- 1-(1-theta.hat)^(data$interval)

  #= Deviance
  Deviance[m] <- -2*sum(dbinom(data$statut_tp1,1,theta.prim,1))
}

ModelComp <- as.data.frame(matrix(NA,nrow=nmodels+1,ncol=2))
names(ModelComp) <- c("Model","Deviance")
ModelComp$Model <-  c(paste("-",varnames,sep=""),"AllVar","Null","Saturated")
ModelComp$Deviance <- c(round(apply(Deviance,1,mean),3),0)

write.table(ModelComp,"./outputs/model/model_comparison.txt",sep="\t",row.names=FALSE,quote=FALSE)

#=====================
#= Variable importance

VarImp <- as.data.frame(matrix(NA,nrow=nmodels-1,ncol=3))
names(VarImp) <- c("Model","Deviance","VarImp")
VarImp$Model <- c(paste("-",varnames,sep=""),"AllVar")
VarImp$Deviance <- round(apply(Deviance,1,mean),3)[1:(nmodels-1)]
VarImp$VarImp <- (VarImp$Deviance-VarImp$Deviance[VarImp$Model=="AllVar"])

write.table(VarImp,"./outputs/model/variable_importance.txt",sep="\t",row.names=FALSE,quote=FALSE)

#=====================
#= Parameter estimates

#= AllVar model
model.spatial <- deforestation(formula=f[[nvar+1]], data=data, interval=data$interval,
                               verbose=1, seed=1234, tune=1, burnin=1000, mcmc=1000, thin=1)

#= Parameter estimates
MCMC <- as.matrix(model.spatial$mcmc)
gamma.hat <- apply(MCMC,2,mean)
sd.hat <- apply(MCMC,2,sd)
q25.hat <- apply(MCMC,2,quantile,0.025)
q75.hat <- apply(MCMC,2,quantile,0.975)
sign0 <- rep(0,length(q75.hat))
sign0[q25.hat*q75.hat>0] <- 1

Estimates <- data.frame(Parameter=colnames(MCMC),gamma=gamma.hat,sd=sd.hat,sign=sign0,row.names=NULL)

write.table(Estimates,"./outputs/model/gamma_full.txt",sep="\t",row.names=FALSE,quote=FALSE)

#=======================================
#= Model with only selected variables

varnames <- c("as.factor(prot_area)","altitude","dist_road","dist_patch",
              "dist_fedge","as.factor(fragindex)","dens")
f <- formula(paste("statut_tp1~",paste(varnames,collapse="+"),sep=""))

model.spatial <- deforestation(formula=f, data=data, interval=data$interval,
                               verbose=1, seed=1234, tune=1, burnin=2000, mcmc=4000, thin=4)

sink("./outputs/model/gamma_selected_coda_output.txt")
print(summary(model.spatial$mcmc))
sink()

pdf("./outputs/model/MCMC_posteriors.pdf")
plot(model.spatial$mcmc)
dev.off()

MCMC <- as.matrix(model.spatial$mcmc)
gamma.hat <- data.frame(gamma=apply(MCMC,2,mean),Parameter=colnames(MCMC),row.names=NULL)

write.table(gamma.hat,"./outputs/model/gamma_selected.txt",sep="\t",row.names=FALSE,quote=FALSE)

#================================
#= Cross-validation and bootstrap

#= Indexes
FOM <- vector()
Sensitivity <- vector()
Specificity <- vector()
OA <- vector() # Overall accuracy
TSS <- vector()
kappa <- vector()

#= Loop on models                 
for (c in 1:ncross) {

  #= subset of data 70/30
  nobs <- nrow(data)
  train.cells <- sort(sample(nobs,size=floor(nobs*70/100),replace=FALSE))
  nobs.train <- length(train.cells)
  test.cells <- c(1:nobs)[-sort(train.cells)]
  nobs.test <- length(test.cells)
  data.train <- data[train.cells,]
  data.test <- data[test.cells,]

  model.spatial <- deforestation(formula=f, data=data.train, interval=data.train$interval,
                                 verbose=1, seed=1234, tune=1, burnin=1000, mcmc=500, thin=1) 
  MCMC <- as.matrix(model.spatial$mcmc)

  #= theta.prim.mean for training data
  data.3 <- rbind(data.1,data.train)
  theta.prim.mean <- sum(data.3$statut_tp1==1)/(sum(data.3$statut_tp1==1)+sum(data.3$statut_tp1==0)) # ML estimates of the deforestation rate
  
  #= Parameter estimates
  gamma.hat <- apply(MCMC,2,mean)

  #= Predictions on test data
  mf.fixed <- model.frame(formula=f,data=data)
  X <- model.matrix(attr(mf.fixed,"terms"),data=mf.fixed)
  theta.hat <- inv.logit(X%*%gamma.hat)
  theta.prim <- 1-(1-theta.hat)^(data$interval)
  theta.prim.test <- theta.prim[test.cells]

  #= Number of pixels to be deforested
  pred <- rep(0,nobs.test)
  ndefor <- round(nobs.test*theta.prim.mean) # number of pixels to be deforested, just for check
  threshold <- quantile(theta.prim.test,1-theta.prim.mean) # probability threshold
  Which <- which(theta.prim.test>=threshold) # pixels with the highest probability of deforestation
  pred[Which] <- 1
  data.test$pred <- pred

  #= Contingency table (pred/obs)
  n00 <- sum(data.test$pred==0&data.test$statut_tp1==0)
  n11 <- sum(data.test$pred==1&data.test$statut_tp1==1)
  n01 <- sum(data.test$pred==0&data.test$statut_tp1==1)
  n10 <- sum(data.test$pred==1&data.test$statut_tp1==0)

  #= Indexes
  OA[c] <- (n11+n00)/(n11+n10+n00+n01)
  FOM[c] <- n11/(n11+n10+n01)
  Sensitivity[c] <- n11/(n11+n01)
  Specificity[c] <- n00/(n00+n10)
  TSS[c] <- Sensitivity[c]+Specificity[c]-1
  kappa[c] <- Kappa(matrix(c(n00,n01,n10,n11),nrow=2))$Unweighted[1]

}

PredPower <- as.data.frame(matrix(NA,nrow=1,ncol=6))
names(PredPower) <- c("OA","FOM","Sen","Spe","TSS","Kappa")

PredPower$OA <- mean(OA)
PredPower$FOM <- mean(FOM)
PredPower$Sen <- mean(Sensitivity)
PredPower$Spe <- mean(Specificity)
PredPower$TSS <- mean(TSS)
PredPower$Kappa <- mean(kappa)

write.table(PredPower,"./outputs/model/model_predictive_power.txt",sep="\t",row.names=FALSE,quote=FALSE)

#===================================================== End of script ===================================================#
