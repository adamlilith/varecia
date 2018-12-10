#===========================================================
# R/GRASS script: Preparing the data-set for modelling deforestation.
#
# Copyright (C) March 2012:
# Clovis Grinand <clovis@goodplanet.org>
# Ghislain Vieilledent <ghislain.vieilledent@cirad.fr>
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
# This R/GRASS script describes the steps to follow in order to obtain
# the data-set for the deforestation model using R and GRASS GIS softwares.
#
# User must open GRASS with location=phcfM_SM and mapset=study_area_4
#==========================================================

#==========================================================
#= The mapsets PERMANENT and study_area_4 already includes the following raw data layers:

#= VECTORS
#-------------------------
# In the PERMANENT mapset
#-------------------------
# instat_fkt@PERMANENT: Fokontany (sub-town level) boundaries and socio-economic data
# instat_frs@PERMANENT: Firaisina (town level) boundaries
# protected_areas_ANGAP@PERMANENT: Protected areas managed by the ANGAP
# main_roads@PERMANENT: Main roads of Madagascar
# main_towns@PERMANENT: Main towns of Madagascar
# Mada@PERMANENT: Madagascar boundaries
# PSs@PERMANENT: Project sites
# SAs@PERMANENT: Study areas
#----------------------------
# In the study_area_4 mapset
#----------------------------
# PS4@study_area_4: Project site #4
# SA4@study_area_4: Study area #4

#= RASTERS
#-------------------------
# In the PERMANENT mapset
#-------------------------
# altitude@PERMANENT: altitude
#----------------------------
# In the study_area_4 mapset
#----------------------------
# fcc@study_area_4: deforestation (forest cover change) on the study area #4
# cloudnshade@study_area_4: clouds and shades on satellites images for the study area #4
# ACD@study_area_4: aboveground carbon density (Mg.ha-1) on the project site #4
# date_t0@study_area_4: dates of the satellite images composing the mosaic at t0 (format: YYYYMMDD)
# date_t1@study_area_4: dates of the satellite images composing the mosaic at t1 (format: YYYYMMDD)
# date_t2@study_area_4: dates of the satellite images composing the mosaic at t2 (format: YYYYMMDD)

#===========================================================
#= In the GRASS console, the user must call R: \"> R\" and can execute the
#= following R/GRASS code

#= Set working directory
# setwd("/path/to/1_deforestation_dataset_SA4_windows.R")

library(spgrass6)

#////////////////////////////////////////////////////////////////////////////////////////////////
#================================================================================================
#
# Step 1: Creating layers of data and explicative variables
#
#================================================================================================
#////////////////////////////////////////////////////////////////////////////////////////////////

library(raster)
library(rgdal)

#= Region settings
shell("g.region vect=SA4 res=30 -pa")

#=======================================================
#= TIME-INTERVAL BETWEEN LAND-COVER OBSERVATIONS

#= Exporting date rasters for use with raster() function in raster package
shell("r.out.gdal input=date_t0 output=./data/date_t0.bil type=UInt32 nodata=9999")
shell("r.out.gdal input=date_t1 output=./data/date_t1.bil type=UInt32 nodata=9999")
shell("r.out.gdal input=date_t2 output=./data/date_t2.bil type=UInt32 nodata=9999")

#= Import in R with raster() function
date.t0 <- raster("./data/date_t0.bil")
date.t1 <- raster("./data/date_t1.bil")
date.t2 <- raster("./data/date_t2.bil")

#= Function to compute time-interval on large rasters
f.difftime <- function (d.start,d.end,ofile) {
  out <- raster(d.start)
  bs <- blockSize(out)
  out <- writeStart(out, ofile, format="EHdr", overwrite=TRUE)
  for (i in 1:bs$n) {
    # a: first date
    a <- getValues(d.start, row=bs$row[i], nrows=bs$nrows[i])
    a1 <- as.character(a)
    a2 <- paste(substr(a1,1,4),substr(a1,5,6),substr(a1,7,8),sep="-")
    a2[a2=="NA-NA-NA"] <- NA
    ua2 <- unique(a2)
    lua2 <- length(unique(a2))
    # b: second date
    b <- getValues(d.end, row=bs$row[i], nrows=bs$nrows[i])
    b1 <- as.character(b)
    b2 <- paste(substr(b1,1,4),substr(b1,5,6),substr(b1,7,8),sep="-")
    b2[b2=="NA-NA-NA"] <- NA
    ub2 <- unique(b2)
    lub2 <- length(unique(b2))
    # difftime
    v <- vector()
    v[is.na(a2) | is.na(b2)] <- NA
    for (j in 1:lua2) {
      for (k in 1:lub2) {
        if (!is.na(ua2[j]) & !is.na(ub2[k])) {
          v[a2==ua2[j] & b2==ub2[k]] <- as.numeric(difftime(ub2[k],ua2[j],units="days")/365.25)
        }
        writeValues(out, v, bs$row[i])
      }
    }
  }
  out <- writeStop(out)
  return(out)
}

# Time-intervals
interval.01 <- f.difftime(date.t0,date.t1,"./data/interval_01.bil") # Take some time to run (~30 secondes)
interval.12 <- f.difftime(date.t1,date.t2,"./data/interval_12.bil") # Take some time to run (~30 secondes)

# Import into GRASS
shell("r.in.gdal input=./data/interval_01.bil output=interval_01 --overwrite")
shell("r.in.gdal input=./data/interval_12.bil output=interval_12 --overwrite")

#=======================================================
#= EXPLICATIVE VARIABLES

#=======================================================
#= Landscape factors

#= Forest and deforestation at t0, t1 and t2
shell("r.mapcalc \"forest_t0=if(fcc==111 || fcc==113 || fcc==133,1,null())\"")
shell("r.mapcalc \"forest_t1=if(fcc==111 || fcc==113,1,null())\"")
shell("r.mapcalc \"forest_t2=if(fcc==111,1,null())\"")
shell("r.mapcalc \"defor_01=if(fcc==133,1,null())\"")
shell("r.mapcalc \"defor_12=if(fcc==113,1,null())\"")

#= Distance to deforested patch between 2000-2005
shell("r.cost -k input=AreaOne output=dist_defor_01 start_rast=defor_01 --overwrite")
shell("r.mapcalc \"dist_defor_01_m=dist_defor_01*(ewres()+nsres())/2.0\"")

#= Distance to forest edge
# t0
shell("r.mapcalc \"nonforest_t0=if(isnull(forest_t0),1,null())\"")
shell("r.cost -k input=AreaOne output=dist_fedge_t0 start_rast=nonforest_t0 --overwrite")
shell("r.mapcalc \"dist_fedge_t0_m=dist_fedge_t0*(ewres()+nsres())/2.0\"")
# t1
shell("r.mapcalc \"nonforest_t1=if(isnull(forest_t1),1,null())\"")
shell("r.cost -k input=AreaOne output=dist_fedge_t1 start_rast=nonforest_t1 --overwrite")
shell("r.mapcalc \"dist_fedge_t1_m=dist_fedge_t1*(ewres()+nsres())/2.0\"")

#= Fragmentation index (with a path to r.forestfrag GRASS script)
shell("sh ./fragindex/r.forestfrag -r input=forest_t0 window=5") # ouput --> fragindex_run raster
shell("g.copy rast=fragindex_run,fragindex_t0 --overwrite")
shell("sh ./fragindex/r.forestfrag -r input=forest_t1 window=5") # ouput --> fragindex_run raster
shell("g.copy rast=fragindex_run,fragindex_t1 --overwrite")

#=======================================================
#= Transport factors

#= Distance to main road
shell("v.to.rast input=main_roads@PERMANENT output=main_roads use=cat --overwrite")
shell("r.cost -k input=AreaOne output=dist_road start_rast=main_roads --overwrite")
shell("r.mapcalc \"dist_road_m=dist_road*(ewres()+nsres())/2.0\"")

#= Distance to main town
shell("v.to.rast input=main_towns@PERMANENT output=main_towns use=cat --overwrite")
shell("r.cost -k input=AreaOne output=dist_town start_rast=main_towns --overwrite")
shell("r.mapcalc \"dist_town_m=dist_town*(ewres()+nsres())/2.0\"")

#=======================================================
#= Socio-economic factors

#= Rasterization of Fokontany socio-economical factors
shell("v.to.rast input=instat_fkt output=fkt_id use=attr column=C_FKT --overwrite")
shell("v.to.rast input=instat_fkt output=fkt_pop use=attr column=F21 --overwrite")
shell("v.to.rast input=instat_fkt output=fkt_date use=attr column=F24 --overwrite")
shell("v.to.rast input=instat_fkt output=fkt_area use=attr column=AIRE_KM2 --overwrite ")
shell("v.to.rast input=instat_fkt output=fkt_mines use=attr column=Mines --overwrite")
shell("v.to.rast input=instat_fkt output=fkt_chemfert use=attr column=ChemFert --overwrite")
shell("v.to.rast input=instat_fkt output=fkt_cattle use=attr column=Cattle --overwrite")
shell("v.to.rast input=instat_fkt output=fkt_poors use=attr column=Poors --overwrite")

#=======================================================
#= Land-policy factors

#= Rasterization of the protected areas
shell("v.to.rast input=protected_areas_ANGAP output=protected_areas_ANGAP use=cat --overwrite")
shell("r.mapcalc \"prot_area=if(isnull(protected_areas_ANGAP),0,1)\"")

#////////////////////////////////////////////////////////////////////////////////////////////////
#================================================================================================
#
# Step 2: Building the data-set
#
#================================================================================================
#////////////////////////////////////////////////////////////////////////////////////////////////

#= Number of random points for each date
nsample <- 20000 # Set a lower number (1000) if you want to test the code quickly

#============================================================
#= Data-set for deforestation process between t0 and t1

#= Data points are sampled outside clouds and shade on satellite images -> MASK
shell("r.mapcalc \"MASK_t=isnull(cloudnshade)\"")
shell("r.mapcalc \"MASK=if(MASK_t==0,null(),1)\"")

#=  Random sample of \"nsample\" points in forest at t0
shell(paste("r.random forest_t0 vect=pts_t0 n=",nsample," --overwrite",sep=""))

#= R vectors with raster and column names 
raster <- c("fcc","interval_01",
            "fkt_id","fkt_date","fkt_pop","fkt_area","fkt_mines","fkt_chemfert","fkt_cattle","fkt_poors",
            "prot_area","altitude","dist_defor_01_m","dist_fedge_t0_m", # Remember that distance to previous deforestation can\"t be used for period 1
            "dist_road_m","dist_town_m","fragindex_t0","date_t0")

colnames <- c("fcc","interval",
             "fkt_id","fkt_date","fkt_pop","fkt_area","fkt_mines","fkt_fert","fkt_cattle","fkt_poors",
             "prot_area","altitude","dist_patch","dist_fedge","dist_road","dist_town","fragindex","date_image")

#= Adding columns to the random points layer
shell("v.db.addcol map=pts_t0 \"columns=fcc INT, interval DOUBLE, fkt_id INT, fkt_date INT, fkt_pop INT, fkt_area DOUBLE, fkt_mines INT, fkt_fert INT, fkt_cattle INT, fkt_poors INT, prot_area INT, altitude INT, dist_patch INT, dist_fedge INT,  dist_road INT, dist_town INT, fragindex INT, date_image INT\"")
	
#= Extracting factor values for each sample point (Take some time to run: ~15 minutes)
for (i in 1:length(colnames)) {
  shell(paste("v.what.rast vector=pts_t0 rast=",  raster[i], " column=", colnames[i], sep=""))
}

#============================================================
#= Data-set for deforestation process between t1 and t2

#=  Random sample of \"nsample\" points in forest at t1
shell(paste("r.random forest_t1 vect=pts_t1 n=",nsample," --overwrite",sep="")) 

#= R vector with raster names 
raster <- c("fcc","interval_12",
            "fkt_id","fkt_date","fkt_pop","fkt_area","fkt_mines","fkt_chemfert","fkt_cattle","fkt_poors",
            "prot_area","altitude","dist_defor_01_m","dist_fedge_t1_m",
            "dist_road_m","dist_town_m","fragindex_t1","date_t1")

#= Adding columns to the random points layer
shell("v.db.addcol map=pts_t1 \"columns=fcc INT, interval DOUBLE, fkt_id INT, fkt_date INT, fkt_pop INT, fkt_area DOUBLE, fkt_mines INT, fkt_fert INT, fkt_cattle INT, fkt_poors INT, prot_area INT, altitude INT, dist_patch INT, dist_fedge INT,  dist_road INT, dist_town INT, fragindex INT, date_image INT\"")
	
#= Extracting factor values for each sample point (Take some time to run: ~15 minutes)
for (i in 1:length(colnames)) {
  shell(paste("v.what.rast vector=pts_t1 rast=",  raster[i], " column=", colnames[i], sep=""))
}

#= Removing mask
shell("r.mask -r")

#============================================================
#= Combining the data-sets for the two periods and exporting

#= Import of the data-sets in R
pts.t0 <- cbind("period"=c(1),readVECT6("pts_t0", ignore.stderr=TRUE, plugin=FALSE)@data)
pts.t1 <- cbind("period"=c(2),readVECT6("pts_t1", ignore.stderr=TRUE, plugin=FALSE)@data)

# Transforming fcc class value in 0 forest and 1 deforested
pts.t0$statut_tp1[pts.t0$fcc==111 || pts.t0$fcc==113] <- 0
pts.t0$statut_tp1[pts.t0$fcc==133] <- 1
pts.t1$statut_tp1[pts.t1$fcc==111] <- 0
pts.t1$statut_tp1[pts.t1$fcc==113] <- 1

# Whole data for the two periods
data.set <- rbind(pts.t0,pts.t1)

# Excluding rows with missing values
na.in.row <- function (x) {sum(is.na(x))>=1}
row.with.na <- apply(data.set,1,na.in.row)
data.set.2 <- data.set[!row.with.na,]

# Exporting the data-set
write.table(data.set.2,file="./data/1_data_deforestation_SA4.txt",sep="\t",quote=FALSE,row.names=FALSE)

#============================================================
#= Cleaning GRASS mapset and output folders

# Uncomment the following lines if you want to clean the mapset and output folder

## #= Rasters to be removed
## lrast <- c("AreaOne","defor_01","defor_12",
##            "dist_cns","dist_cns_m",
##            "dist_defor_01","dist_defor_01_m","dist_fedge_t0",         
##            "dist_fedge_t0_m","dist_fedge_t1","dist_fedge_t1_m",     
##            "dist_road","dist_road_m","dist_town",             
##            "dist_town_m","fkt_area","fkt_cattle",            
##            "fkt_chemfert","fkt_date","fkt_id",                
##            "fkt_mines","fkt_poors","fkt_pop",               
##            "forest_t0","forest_t1","forest_t1",             
##            "forest_t2","forest_t2","fragindex_run",         
##            "fragindex_t0","fragindex_t1","interval_01",           
##            "interval12","interval_12","main_roads",            
##            "main_towns","nonforest_t0","nonforest_t1",          
##            "prot_area","protected_areas_ANGAP","PS4",                   
##            "SA4")

## for (i in 1:length(lrast)) {
##   shell(paste("g.remove rast=",lrast[i],sep=""))
## }

## #= Vectors
## lvect <- c("pts_t0,pts_t1")

## for (i in 1:length(lvect)) {
##   shell(paste("g.remove vect=",lvect[i],sep=""))
## }

## #= Output folder
## shell("rm ./data/1_data_deforestation_SA4.txt")
## shell("rm ./data/data_*")
## shell("rm ./data/interval_*")

#= End of script
#=============================================================================================#








