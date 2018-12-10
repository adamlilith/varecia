#==========================================================
# R/GRASS script: Forecasting deforestation.
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
# This R script describes the steps to follow in order to forecast the
# deforestation model using R and GRASS GIS softwares.
#
# User must open GRASS with location=phcfM_SM and mapset=study_area_4
#==========================================================

#= Set working directory
# setwd("/path/to/3_deforestation_forecast_PS4_windows.R")

#= Loading libraries
library(spgrass6)

#= inv.logit function
inv.logit <- function(x, min=0, max=1) {
  p <- exp(x)/(1+exp(x))
  p <- ifelse(is.na(p) & !is.na(x), 1, p) # fix problems with +Inf
  return(p*(max-min)+min)
}

#= Remove potential mask
## shell("r.mask -r")

#= Population effect
pop <- TRUE # If you don't want the population effect to be taken into account, set pop to FALSE

#= Create \"forecast\" and \"temp\" folder
shell("rm -r ./outputs/forecast")
shell("rm -r ./outputs/temp")
shell("mkdir .\\outputs\\forecast")
shell("mkdir .\\outputs\\temp")

# Rasterize SA4 and PS4
shell("v.to.rast --o input=SA4 output=SA4 column=cat")
shell("v.to.rast --o input=PS4 output=PS4 column=cat")

#= Some statistics
shell("g.region vect=SA4 res=30 -ap")
shell("r.mapcalc \"forest2010=if(fcc==111,1,null())\"")
shell("r.report SA4 units=h > ./outputs/forecast/area_sa.txt")
shell("r.report PS4 units=h > ./outputs/forecast/area_ps.txt")
shell("r.mask input=SA4")
shell("r.report forest2010 units=h > ./outputs/forecast/forest_sa.txt")
shell("r.mask -o input=PS4")
shell("r.report forest2010 units=h > ./outputs/forecast/forest_ps.txt")
shell("r.mask -o input=forest2010")
shell("r.univar ACD > ./outputs/forecast/acd_ps.txt")
shell("r.mask -r")

#= Computing the probability of deforestation for each pixel
#
# model is theta_i=inv.logit(gamma_0+gamma_1j*fragindex_it+gamma_2*aire_prot_it+gamma_3*dens_it)
# with inv.logit=exp(x)/(1+exp(x))=1/(1+exp(-x))

#= Parameter values for probability of deforestation
gamma0 <- -1.62592922974288	# (Intercept)
gamma1 <- -0.349459614738086	# as.factor(prot_area)1: presence
gamma2 <- -0.000956401028177047	# altitude
gamma3 <- -5.31044050617831e-06	# dist_road
gamma4 <- -0.0027068646371299	# dist_patch
gamma5 <- -0.00398188606513584	# dist_fedge
gamma61 <- 0                    # as.factor(fragindex)1: patch
gamma62 <- -0.0495266548790051	# as.factor(fragindex)2: transitional
gamma63 <- -0.419137256802682	# as.factor(fragindex)3: perforated
gamma64 <- -0.307849906952335	# as.factor(fragindex)4: edge
gamma65 <- -0.471731607426946	# as.factor(fragindex)5: interior
gamma7 <- 0.00620052366273907	# dens

#= Population growth rate
rho <- 0.03389765 # population growth rate

#= Preparing rasters of covariates at time=t
year <- c(2011:2030)
nyear <- length(year)

#= Parameter value for intensity of deforestation
beta <- read.table("./outputs/model/beta.txt",sep="\t",header=TRUE)
beta0 <- beta$beta[1]
beta1 <- beta$beta[2]

#= Project region
shell("g.region -ap vect=PS4 res=30")
shell("r.mask input=PS4")
shell("r.mapcalc \"forest_2010=forest2010\"")

# Initial forest raster and covariates
shell("r.mapcalc \"forest_run=forest_2010\"") # Put here the starting date for forest
shell("r.mapcalc \"defor_patch_run=if(fcc==113,1,null())\"")  # Patches of deforestation between t1 and t2

# Raster of K (population at year=1993)
shell(paste("r.mapcalc \"K_pop=fkt_pop/exp(",rho,"*(fkt_date-1993))\"",sep=""))

# AreaOne for r.cost and distance computation
shell("r.mapcalc \"AreaOne=1\"")

#============================
# Loop on years of simulation
for (i in 1:nyear) {

  #= Computing fragindex (time consuming)
  shell("sh ./fragindex/r.forestfrag -r input=forest_run window=5") # ouput --> fragindex_run raster

  #= Computing distance to previously deforested patch
  shell("r.cost -k input=AreaOne output=dist_defor_patch start_rast=defor_patch_run --overwrite")
  shell("r.mapcalc \"dist_defor_patch_run=dist_defor_patch*(ewres()+nsres())/2\"")
  
  #= Computing distance to forest edge
  shell("r.mapcalc \"nonforest_run=if(isnull(forest_run),1,null())\"")
  shell("r.cost -k input=AreaOne output=dist_fedge start_rast=nonforest_run --overwrite")
  shell("r.mapcalc \"dist_fedge_run=dist_fedge*(ewres()+nsres())/2\"")

  #= Decomposing fragindex in 5 rasters (patch, transitional, edge, perforated, interior)
  shell("r.mapcalc \"patch_run=if(fragindex_run==1,1,0)\"")
  shell("r.mapcalc \"transitional_run=if(fragindex_run==2,1,0)\"")
  shell("r.mapcalc \"edge_run=if(fragindex_run==3,1,0)\"")
  shell("r.mapcalc \"perforated_run=if(fragindex_run==4,1,0)\"")
  shell("r.mapcalc \"interior_run=if(fragindex_run==5,1,0)\"")

  #= Computing density on each pixel and average density on PS
  # On each pixel
  if (pop) {
    shell(paste("r.mapcalc \"pop_run=K_pop*exp(",rho,"*(",year[i]-1,"-1993))\"",sep=""))
  }
  if (!pop) {
    shell(paste("r.mapcalc \"pop_run=K_pop*exp(",rho,"*(",2010,"-1993))\"",sep=""))
  }
  shell("r.mapcalc \"dens_run=pop_run/fkt_area\"")
  # Average on PS
  shell("r.average base=PS4 cover=dens_run out=average_dens_run --o")
  shell("r.univar average_dens_run > ./outputs/temp/average_dens_run.txt")
  S <- scan(file="./outputs/temp/average_dens_run.txt",what="character")
  dens.mean.run <- as.numeric(S[25])
  
  #= Intensity of deforestation
  theta.intensity.run <- inv.logit(beta0+beta1*dens.mean.run)

  #= Proba of deforestation
  shell(paste("r.mapcalc \"theta_run=if(forest_run==1,1/(1+exp(-(",gamma0,"+"
               ,gamma1,"*prot_area+"
               ,gamma2,"*altitude+"
               ,gamma3,"*dist_road_m+"
               ,gamma4,"*dist_defor_patch_run+"
               ,gamma5,"*dist_fedge_run+"
               ,gamma61,"*patch_run+"
               ,gamma62,"*transitional_run+"
               ,gamma63,"*edge_run+"
               ,gamma64,"*perforated_run+"
               ,gamma65,"*interior_run+"
               ,gamma7,"*dens_run))),null())\"",sep=""))
  
  #= Threshold using quantiles
  shell(paste("r.quantile input=theta_run quantiles=1 percentiles=",100*(1-theta.intensity.run)," > ./outputs/temp/threshold_run.txt",sep=""))
  S <- scan("./outputs/temp/threshold_run.txt",what=double(),sep=":")
  threshold.run <- S[3]

  #= Deforesting following probability threshold
  shell(paste("r.mapcalc \"defor_patch_run=if(theta_run>=",threshold.run,",1,null())\"",sep=""))
  shell("r.mapcalc \"forest_run_t=forest_run\"")
  shell(paste("r.mapcalc \"forest_run=if(theta_run>=",threshold.run,",null(),forest_run_t)\"",sep=""))
  
  #= Evolution of the variables with time
  if (i==1) {
    cat(paste(theta.intensity.run,"\n",sep=""), file="./outputs/forecast/evol_theta_intensity.txt", append=FALSE)
    cat(paste(dens.mean.run,"\n",sep=""), file="./outputs/forecast/evol_dens_mean.txt", append=FALSE)
  }
  if (i>1) {
    cat(paste(theta.intensity.run,"\n",sep=""), file="./outputs/forecast/evol_theta_intensity.txt", append=TRUE)
    cat(paste(dens.mean.run,"\n",sep=""), file="./outputs/forecast/evol_dens_mean.txt", append=TRUE)
  }
  shell("r.mask -o input=forest_run")
  shell("r.univar altitude >> ./outputs/forecast/evol_time_altitude.txt")
  shell("r.univar fragindex_run >> ./outputs/forecast/evol_time_fragindex.txt")
  shell("r.univar ACD >> ./outputs/forecast/evol_time_ACD.txt")
  shell("r.mask -o input=PS4")

  #= Output
  shell(paste("r.mapcalc \"forest_",year[i],"=forest_run\"",sep=""))
  
} # End for loop on years of simulation

#===============================
# Simplifying variable evolution files with awk
shell("gawk '$1 == \"mean:\" {print}' ./outputs/forecast/evol_time_altitude.txt > ./outputs/forecast/evol_altitude.txt")
shell("gawk '$1 == \"mean:\" {print}' ./outputs/forecast/evol_time_fragindex.txt > ./outputs/forecast/evol_fragindex.txt")
shell("gawk '$1 == \"mean:\" {print}' ./outputs/forecast/evol_time_ACD.txt > ./outputs/forecast/evol_ACD.txt")
shell("rm ./outputs/forecast/evol_time*") # Remove temporary 'evolution' files
shell("rm ./outputs/temp/*_run.txt") # Remove temporary 'run' files

#===============================
# Changing the color for rasters
for (i in 1:nyear) {
  shell(paste("echo 1 26:152:80 | r.colors forest_",year[i]," color=rules",sep=""))
}

#================================================
# Quantifying carbon released from deforestation
# On the PHCF project zone

#= Raster of Carbon in Mg (ACD in Mg.ha-1)
shell("r.mapcalc \"Carbon=ACD*(ewres()*nsres())/10000\"") # 1 pixel is 30*30=900 m^2

#= Released carbon
shell("r.mask -o input=forest_2010")
shell("r.sum rast=Carbon >> ./outputs/forecast/C_evolution.txt")
for (i in 1:nyear) {
  shell(paste("r.mask -o input=forest_",year[i],sep=""))
  shell("r.sum rast=Carbon >> ./outputs/forecast/C_evolution.txt")
}
#= Removing MASK
shell("r.mask -r")

#= Carbon evolution graphics

#= Importing C_evolution into R
Data.C <- read.table(file="./outputs/forecast/C_evolution.txt",header=FALSE,sep="=")
names(Data.C) <- c("Year","C_stock")
Data.C$Year <- c(2010,year)
Data.C$C_release <- 0
for (i in 2:(nyear+1)) {
  Data.C$C_release[i] <- -(Data.C$C_stock[i]-Data.C$C_stock[i-1])
  Data.C$C_release[i] <- Data.C$C_release[i]+Data.C$C_release[i-1]
}
Data.C$CO2_release <- Data.C$C_release*44/12 # 1 metric ton of carbon dioxide equivalent (tCO2e) = 44/12 metric tons carbon

#= Loop on years
for (i in 1:nrow(Data.C)) {
  png(paste(file="./outputs/forecast/CO2_released_",Data.C$Year[i],".png",sep=""),width=500,height=500)
  par(cex=1.4,mar=c(5,6,3,1))
  plot(c(1:i),Data.C$CO2_release[c(1:i)],
       type="l",
       main=paste("Year: ",Data.C$Year[i],sep=""),
       xlim=c(1,nrow(Data.C)),
       ylim=c(0,4000000),
       xlab="Year",
       ylab="CO2 released \nfrom deforestation (GT)",
       cex.lab=1.4,
       cex.main=1.4,
       axes=FALSE)
  axis(1,at=seq(1,nrow(Data.C),1),labels=rep("",nrow(Data.C)),tcl=0.5,line=0.25)
  axis(1,at=seq(1,nrow(Data.C),5),labels=seq(2010,2030,by=5),cex.axis=1.4,tcl=0.5,line=0.25)
  axis(2,at=seq(from=0,to=4000000,by=1000000),labels=seq(from=0,to=4,by=1),cex.axis=1.4,tcl=0.5,line=0.25)
  if (i==1) {
    points(1,Data.C$CO2_release[1],pch=16,cex=1.5,col="black")
  }
  if (i>1 & i<=11) {
    points(1,Data.C$CO2_release[1],pch=16,cex=1.5,col="black")
    points(c(2:i),Data.C$CO2_release[c(2:i)],pch=16,cex=1.5,col="orange")
  }
  if (i>11 & i<=21) {
    points(1,Data.C$CO2_release[1],pch=16,cex=1.5,col="black")
    points(c(2:i),Data.C$CO2_release[c(2:i)],pch=16,cex=1.5,col="orange")
    points(c(12:i),Data.C$CO2_release[c(12:i)],pch=16,cex=1.5,col="red")
  }
  text(x=11,y=4000000,label=paste(round(Data.C$CO2_release[i])," T",sep=""),cex=2)
  if (Data.C$Year[i]>=2020) {
    text(x=11,y=Data.C$CO2_release[Data.C$Year==2020]+1250000,
         label=paste(round(Data.C$CO2_release[Data.C$Year==2020])," T",sep=""))
    segments(x0=11,y0=Data.C$CO2_release[Data.C$Year==2020]+1100000,
             x1=11,y1=Data.C$CO2_release[Data.C$Year==2020]+150000)
  }
  dev.off()
}

#===================================================
# Animation on the whole project area

#= Set region to project area
shell("g.region -ap vect=PS4 res=30")

#= Marking deforestation with red pixels (value:0)

# Difference in land cover compared to the start
shell("r.mapcalc \"deforest_PS_2010=forest_2010\"")
# First 10 years
for (i in 2:11) {
  shell(paste("r.mapcalc \"deforest_PS_",Data.C$Year[i],"=if(forest_",Data.C$Year[1],"==1 && isnull(forest_",Data.C$Year[i],"),2,forest_",Data.C$Year[i],")\"",sep=""))
}
# Second 10 years
shell("r.mapcalc \"forest_2020_t=forest_2020\"")
shell("r.mapcalc \"forest_2020=if(isnull(forest_2020_t),0,forest_2020_t)\"") # This is needed to use if in mapcalc
for (i in 12:21) {
  shell(paste("r.mapcalc \"deforest_PS_",Data.C$Year[i],"=if(forest_",Data.C$Year[1],"==1 && isnull(forest_",Data.C$Year[i],"),2,forest_",Data.C$Year[i],")\"",sep=""))
  shell(paste("r.mapcalc \"deforest_PS_",Data.C$Year[i],"_t=deforest_PS_",Data.C$Year[i],"\"",sep=""))
  shell(paste("r.mapcalc \"deforest_PS_",Data.C$Year[i],"=if(forest_",Data.C$Year[11],"==1 && isnull(forest_",Data.C$Year[i],"),3,deforest_PS_",Data.C$Year[i],"_t)\"",sep=""))
}
shell("r.mapcalc \"forest_2020_t=forest_2020\"")
shell("r.mapcalc \"forest_2020=if(forest_2020_t==0,null(),forest_2020_t)\"") # 0 replaced by null()

# Change color legend
for (i in 1:nrow(Data.C)) {
  shell(paste("(echo 1 26:152:80 && echo 2 orange && echo 3 red) | r.colors deforest_PS_",Data.C$Year[i]," color=rules",sep=""))
}
  
#= Exporting rasters as images
for (i in 1:nrow(Data.C)) {
  shell(paste("r.out.png input=deforest_PS_",Data.C$Year[i]," output=./outputs/forecast/deforest_PS_",Data.C$Year[i],".png",sep=""))
}

#= Concatenate deforestation images with Carbon evolution images
#  (Image Magick required: www.imagemagick.org)
for (i in 1:nrow(Data.C)) {
  shell(paste("montage ./outputs/forecast/deforest_PS_",Data.C$Year[i],".png ./outputs/forecast/CO2_released_",Data.C$Year[i],".png -geometry 500x500+1+1 ./outputs/forecast/defor_PS_CO2_",Data.C$Year[i],".png",sep=""))
}

#= Animated .GIF
shell("mogrify -format gif ./outputs/forecast/defor_PS_CO2_*.png") # Convert images to .gif
shell("convert -delay 100 -loop 1 ./outputs/forecast/defor_PS_CO2_*.gif ./outputs/forecast/film_defor_PS.gif")
shell("mv ./outputs/forecast/defor_PS_CO2_2030.png ./outputs/forecast/PS_CO2_2030.png") # Last image
shell("rm ./outputs/forecast/defor*") # Removing "defor" files
shell("rm ./outputs/forecast/CO2*") # Removing "CO2" files

## #= Uncomment the following lines if you want to clean the mapset and output folder
## system ("g.mremove -fb rast=*_run")
## system ("g.mremove -fb rast=forest_*")
## system ("g.mremove -fb rast=deforest_*")
## system ("g.remove -f rast=AreaOne,Carbon,dist_defor_patch,dist_fedge,K_pop")

#= Removing temp folder
shell("rm -R ./outputs/temp") # Be careful with rm -R command, can delete all your data !!!


#= End of script
#=============================================================================================#
