### Modeling range of Varecia, the Black-and-White Lemur (with Toni Lyn Morelli)
### Adam B. Smith | Missouri Botanical Garden | adamDOTsmithATmobotDOTorg | 2018-02

### This code is version 6 of attempts to model the exposure of two lemurs, Varecia variegata and Varecia rubra, to anticipated climate change and deforestation in Madagascar. The analysis depends on two models, one of deforestation and an ecological niche model using forest cover and climate as predictors.

# source('C:/Ecology/Drive/Research/Varecia (Toni Lyn Morelli)/Versions 06/Code/Varecia Ver 06 Ecological Niche Modeling.r')
# source('H:/Global Change Program/Research/Varecia (Toni Lyn Morelli)/Versions 06/Code/Varecia Ver 06 Ecological Niche Modeling.r')
# source('E:/ecology/Drive/Research/Varecia (Toni Lyn Morelli)/Versions 06/Code/Varecia Ver 06 Ecological Niche Modeling.r')

	setwd('C:/ecology/Drive/Research/Varecia (Toni Lyn Morelli)/Versions 06')
	# setwd('H:/Global Change Program/Research/Varecia (Toni Lyn Morelli)/Versions 06')
	# setwd('E:/ecology/Drive/Research/Varecia (Toni Lyn Morelli)/Versions 06')

################
### CONTENTS ###
################	

### SETUP ###

	### libraries, functions, and definitions ###
	### create masks of Madagascar (entire island) ###
	### collate current environmental data for ecological niche models ###
	
	### collate species data ###
	### create display of occurrence records ###
	### create display of Madagascar with PAs ###
	### match occurrence records with environmental data ###
	### select background sites ###
	### calculate spatial autocorrelation between survey sites ###
	### calculate survey site weights to correct for sampling bias ###
	### collate occurrence records and background sites ###
	### train ecological niche models ###
	
	### select and collate future climate data for projecting ecological niche model ###
	
	### visualize response curves of ecological niche models ###
	### evaluate ecological niche models ###
	
	### predict to current and future conditions ###
	### create ensemble rasters for future climate ecological niche model projections ###
	### calculate mean environmental suitability across entire region and elevation bands ###
	### report values of mean environmental suitability across entire region and elevation bands ###
	### compare elevational distribution of forest and occurrences ###
	### evaluate changes in suitability in protected areas ###
	
	### create fine-scale hillshade raster ###
	### create displays of ecological niche model predictions ###
	### create 3D displays of ecological niche model predictions ###
	### create 3D displays of climate change ###
	### create display range maps by surveyor ###
	
	### niche overlap analysis ###
	

#############################################
### libraries, functions, and definitions ###
#############################################

	memory.limit(memory.limit() * 2^30)
	rm(list=ls())
	gc()

	library(omnibus)
	library(enmSdm)
	library(statisfactory)
	library(legendary)
	library(fasterRaster)
	
	library(sp)
	library(graticule)
	library(raster)
	library(parallel)
	library(dismo)
	library(rgeos)
	library(geosphere)
	
	library(brglm2)
	# library(phcfM)
	
	library(fpCompare)
	library(scales)
	library(rgl)
	# library(rayshader)
	library(tictoc)
	
	say(date())

	###############
	### options ###
	###############
	
		options(stringsAsFactors=FALSE)
		rasterOptions(format='GTiff', overwrite=TRUE)

	#################
	### variables ###
	#################
		
		taxa <- c('genus', 'variegata', 'rubra')
		
		madEaProj <- '+init=epsg:32738'
		
		bioclims <- c(4, 10, 15, 16) # BIOCLIM variables to use

		predictors <- c(paste0('bio', prefix(bioclims, 2)), 'forestFragClass_utm38s') # predictors in ENMs

		longLat <- c('longWgs84', 'latWgs84') # BIOCLIM variables to use
		
		grassDir <- c('C:/OSGeo4W64/', 'grass-7.4.1', 'osgeo4W')

		periods <- c('2041 to 2060', '2061 to 2080') # future time periods corresponding to WORLDCLIM data
		rcps <- c('2pt6', '4pt5', '6pt0', '8pt5') # RCPs for future climate periods
		# rcps <- c('2pt6', '6pt0') # RCPs for future climate periods
		# rcps <- c('4pt5', '8pt5') # RCPs for future climate periods

		### extents of insets for plots of maps
		madFocus1 <- extent(929471, 1043462, 8234420, 8358325)
		madFocus2 <- extent(850269, 944754, 7949353, 8048591)
		madFocus3 <- extent(687756, 751211.0, 7433385, 7553616)

		madFocus1 <- as(madFocus1, 'SpatialPolygons')
		madFocus2 <- as(madFocus2, 'SpatialPolygons')
		madFocus3 <- as(madFocus3, 'SpatialPolygons')

		projection(madFocus1) <- projection(madEaProj)
		projection(madFocus2) <- projection(madEaProj)
		projection(madFocus3) <- projection(madEaProj)

	#################
	### functions ###
	#################

		# convert NAs to zeros
		naToZeroFx <- function(x) ifelse(is.na(x), 0, x)

		# equation for annualized forest loss rate from Puyravaud, J-P. 2003 Standardizing the calculation of the annual rate of deforestation. Forest Ecology and Management 177:593-596.
		# t1, t1	Time 1 and 2
		# A1, A2	Area of forest in times 1 and 2
		compoundInterestLaw <- function(t1, t2, A1, A2) (1 / (t2 - t1)) * log(A2 / A1)
			
		# predict ecological niche model, convert to integer to reduce file size
		predictVareciaEnm <- function(model, predStack, cores=6) {

			# model 	model object
			# preds 	predictor stack
			# cores		nummber of cores

			beginCluster(cores)
				prediction <- clusterR(preds, predict, args=list(model=model, type='response'))
				prediction <- 100 * prediction
				prediction <- clusterR(prediction, round)
			endCluster()
			
			prediction
			
		}
		
		# calculate ensemble mean (across GCMs, say) for projection rasters from ecological niche model, convert to integer to reduce file size
		ensembleVareciaEnmRasters <- function(predictions, cores=6) {
			
			# predictions	raster stack of prediction rasters
			# cores			number of cores
		
			rangeFx <- function(x) max(x) - min(x)
			beginCluster(cores)
				meanPrediction <- clusterR(predictions, calc, args=list(fun=mean))
			endCluster()
			
			meanPrediction
			
		}
			
		# return nice names for predictors
		predictorNice <- function(p, incUnits=FALSE) {
			
			out <- if (p == 'bio04') {
				'Temperature Seasonality'
			} else if (p == 'bio10') {
				'Warmest Quarter Temperature'
			} else if (p == 'bio15') {
				'Precipitation Seasonality'
			} else if (p == 'bio16') {
				'Wettest Quarter Precipitation'
			} else if (p == 'forestFragmentationClass') {
				'Forest Fragmentation Class'
			}
			
			if (incUnits) {
				if (p == 'bio10') {
					out <- paste0(out, ' (deg C)')
				} else if (p == 'bio16') {
					out <- paste0(out, ' (mm)')
				}
			}
				
			out
			
		}
			

		# balance sum of weights of presences and background sites
		balanceWeights <- function(x) {
		
			# df	data frame with at least these fields: 'presBg', 'weight'

			weight <- x$weight
			
			presWeight <- sum(weight[x$presBg == 1])
			bgWeight <- sum(weight[x$presBg == 0])
			
			if (presWeight > bgWeight) {
				weight[x$presBg == 0] <- weight[x$presBg == 0] * (presWeight / bgWeight)
			} else {
				weight[x$presBg == 1]<- weight[x$presBg == 1] * (bgWeight / presWeight)
			}
			
			weight
		
		}

		### convert singlepart spatial polygons to multipart
		multiPart <- function(shape) {

			#vector of number of elementary polygons in each multipart polygon (for calculation of plot order, and initialization of islands)
			nbpoly <- sapply(shape@polygons, function(x){ length(x@Polygons) })

			#initialize the list of single Polygons to be built
			islands <- vector('list', sum(nbpoly))
			
			#initialize the vector plot order of single polygons to be built
			plotorder <- vector('integer',sum(nbpoly))

			#index of current single Polygons
			k <- 0

			#loop on all Polygons of the object

			for (i in 1:length(shape)) {

				#current multiple polygon
				pols <- shape@polygons[[i]]@Polygons
				ID <- shape@polygons[[i]]@ID

				#number of polygons to plot before the current multiple one
				prev <- sum(nbpoly[shape@plotOrder < shape@plotOrder[i]])

				#loop on each elementary polygon of the current multiple polygon

				for (j in 1:length(pols)) {
					k <- k+1
					IDs <- ifelse(length(pols)>1,paste(ID,'-',j,sep=''),ID)
					islands[[k]] <- Polygons(list(pols[[j]]), IDs)
					plotorder[k] <- rank(shape@polygons[[i]]@plotOrder,ties.method="first")[j]+prev
				}
			
			}

			multitoone <- SpatialPolygons(islands,pO=plotorder)
			if (!is.na(proj4string(shape))) proj4string(multitoone) <- proj4string(shape)

			multitoone
			
		}

		### stack *current* predictors that use WGS84 30 arcsec resolution for Madagascar
		#################################################################################
		
		stackMadCurrent_wgs84 <- function() {
		
			yrs <- c(2000, 2005, 2010, 2014)
		
			out <- raster::stack(c(
				'./Data/Topography - WORLDCLIM Ver 1pt4 Rel 3/elevation_wgs84.tif',
				paste0('./Data/Forest - Vieilledent et al 2018/forest', yrs, '_wgs84.tif'),
				paste0('./Data/Climate - WORLDCLIM Ver 2 Rel 1 Jun 2016 (1970-2000)/bio', prefix(bioclims, 2), '.tif')
			))
			
			out
		
		}

		### stack *future* climate predictors that use WGS84 30 arcsec resolution for Madagascar
		########################################################################################

		stackMadFuture_wgs84 <- function(gcm, period, rcp) {
		
			# gcm
			# period
			# rcp
		
			out <- raster::stack(
				paste0('./Data/Climate - WORLDCLIM Ver 2 Rel 1 Jun 2016 (', period, ')/', gcm, ' RCP', rcp, '/bio', prefix(bioclims, 2), '.tif')
			)
			
			out
		
		}

		### stack *current* predictors that use WGS84 30 arcsec resolution for Madagascar
		#################################################################################

		stackMadCurrent_utm38s <- function() {
		
			yrs <- c(2000, 2005, 2010, 2014)
		
			out <- raster::stack(c(
				paste0('./Data/Forest - Vieilledent et al 2018/forest', yrs, '.tif'),
				paste0('./Data/Forest - Vieilledent et al 2018/forestFragClass', yrs, '_utm38s.tif'),
				paste0('./Data/Forest - Vieilledent et al 2018/forestFragConnect', yrs, '_utm38s.tif'),
				paste0('./Data/Forest - Vieilledent et al 2018/forestFragDensity', yrs, '_utm38s.tif')
			))
			
			names(out)[names(out) == 'forest2000'] <- 'forest2000_utm38s'
			names(out)[names(out) == 'forest2005'] <- 'forest2005_utm38s'
			names(out)[names(out) == 'forest2010'] <- 'forest2010_utm38s'
			names(out)[names(out) == 'forest2014'] <- 'forest2014_utm38s'
	
			out
		
		}
		
# say('##################################################')
# say('### create masks of Madagascar (entire island) ###')
# say('##################################################')

	# elev <- raster::getData(name='alt', country='MDG', path='C:/ecology/!Scratch')
	# mask_wgs84 <- elev * 0 + 1
	# names(mask_wgs84) <- 'mask_wgs84'
	# writeRaster(mask_wgs84, './Study Region & Masks/WGS84 30-arcsec Resolution/mask_wgs84', datatype='INT1U')
		
# say('######################################################################')
# say('### collate current environmental data for ecological niche models ###')
# say('######################################################################')		
	
	# mask_wgs84 <- raster('./Study Region & Masks/WGS84 30-arcsec Resolution/mask_wgs84.tif')
	# humidForestBufferMask_utm38s <- raster('./Study Region & Masks/UTM 38S 30-m Resolution/humidForestBufferMask_utm38s.tif')
	
	# say('elevation...')
	
		# dirCreate('./Data/Topography - WORLDCLIM Ver 1pt4 Rel 3')
		# elev <- raster::getData(name='alt', country='MDG', path='C:/ecology/!Scratch')
		# names(elev) <- 'elevation_wgs84'
		# writeRaster(elev, './Data/Topography - WORLDCLIM Ver 1pt4 Rel 3/elevation_wgs84', datatype='INT4S')

	# say('current forest cover...')
	
		# forests <- c('forest2000', 'forest2005', 'forest2010', 'forest2014')

		# for (thisForest in forests) {
			# assign(thisForest, raster(paste0('./Data/Forest - Vieilledent et al 2018/', thisForest, '.tif')))
		# }
		
		# for (thisForest in forests) {
		
			# say(thisForest)
		
			# x <- get(thisForest)
			
			# beginCluster(4)
			# x0 <- clusterR(x, calc, args=list(fun=naToZeroFx))
			# endCluster()
			
			# x0 <- humidForestBufferMask_utm38s * x0
			
			# gc()
			# xFract <- fasterFocal(x0, w=3, fun=mean, na.rm=TRUE, pad=TRUE, cores=3)
			# gc()
		
			# beginCluster(4)
			# xProj <- projectRaster(xFract, mask_wgs84, method='bilinear')
			# endCluster()
			
			# xProj <- mask_wgs84 * xProj
			
			# names(xProj) <- paste0(thisForest, '_wgs84')
			# writeRaster(xProj, paste0('./Data/Forest - Vieilledent et al 2018/', thisForest, '_wgs84'))
			
		# }
		
		
	# say('current climate...')
	
		# dirCreate('./Data/Climate - WORLDCLIM Ver 2 Rel 1 Jun 2016 (1970-2000)')
	
		# for (thisBio in c(bioclims, 1, 12)) {
		
			# rast <- raster(paste0('D:/ecology/Climate/WORLDCLIM Ver 2 Rel June 1 2016/30 arcsec 1970 to 2000/wc2.0_bio_30s_', prefix(thisBio, 2), '.tif'))
			
			# rast <- crop(rast, mask_wgs84)
			# names(rast) <- paste0('bio', prefix(thisBio, 2))
			# writeRaster(rast, paste0('./Data/Climate - WORLDCLIM Ver 2 Rel 1 Jun 2016 (1970-2000)/bio', prefix(thisBio, 2)))
	
		# }
	
# say('############################')
# say('### collate species data ###')
# say('############################')	

	# varecia <- read.csv('./Data/Varecia/02 Master__Varecia_Point_Data_2017 - Removed and Corrected Badly Georeferenced Sites and Zahamena site by EEL on 3-10-2000.csv')

	# say('Adding data from Ratsimbazafy...')
		
		# # points already in WGS84
		# ratWgs84 <- read.csv('./Data/Varecia/Ratsimbazafy Sites in WGS84.csv')
		# varecia <- rbind(varecia, ratWgs84)
	
		# # points in UTM 38S
		# ratUtm38s <- read.csv('./Data/Varecia/Ratsimbazafy Sites in UTM 38S.csv')		
		# ratUtm38sPoints <- SpatialPoints(cbind(ratUtm38s$LONG_UTM38S, ratUtm38s$LAT_UTM38S), proj4=CRS(madEaProj))
		# ratUtm38sPoints <- sp::spTransform(ratUtm38sPoints, getCRS('wgs84', TRUE))
		
		# coords <- coordinates(ratUtm38sPoints)
		# coords <- coords[ , 2:1]
		# coords <- as.data.frame(coords)
		# names(coords) <- c('LAT_WGS84', 'LONG_WGS84')
		
		# ratUtm38s <- cbind(ratUtm38s[ , c('Species', 'Status')], coords, ratUtm38s[ , c('Provider', 'Site', 'Date')])
		# varecia <- rbind(varecia, ratUtm38s)

	# say('Renaming columns...')	
		
		# names(varecia)[names(varecia) %in% 'Species'] <- 'species'
		# names(varecia)[names(varecia) %in% 'Status'] <- 'presAbs'
		# names(varecia)[names(varecia) %in% 'LAT_WGS84'] <- 'latWgs84'
		# names(varecia)[names(varecia) %in% 'LONG_WGS84'] <- 'longWgs84'
		# names(varecia)[names(varecia) %in% 'Provider'] <- 'provider'
		# names(varecia)[names(varecia) %in% 'Site'] <- 'site'
		# names(varecia)[names(varecia) %in% 'Date'] <- 'origDate'
		
	# say('Parsing years...')	

		# varecia$year <- omnibus::yearFromDate(varecia$origDate)
		
		# for (i in 1:nrow(varecia)) {
			# if (!is.na(varecia$year[i])) {
				# if (varecia$year[i] >= 9900 & varecia$year[i] <= 9916) {
					# varecia$year[i] <- varecia$year[i] - 9900 + 2000
				# } else if (varecia$year[i] >= 9997) {
					# varecia$year[i] <- varecia$year[i] - 9900 + 1900
				# }
			# }
		# }
		
		# varecia$year[varecia$origYear == 'after 1989'] <- 1990
	
	# varecia <- varecia[varecia$longWgs84 > 44, ]
	
	# write.csv(varecia, './Data/Varecia/03 Added Ratsimbazafy and Parsed Year.csv', row.names=FALSE)

# say('############################################')
# say('### create display of occurrence records ###')
# say('############################################')	

	# outDir <- './Figures & Tables/Maps of Occurrences/'
	# dirCreate(outDir)

	# # hillshading
	# elev <- raster('./Data/Topography - WORLDCLIM Ver 1pt4 Rel 3/elevation_wgs84.tif')
	# slope <- terrain(elev, 'slope')
	# aspect <- terrain(elev, 'aspect')
	
	# hs <- hillShade(slope, aspect, angle=315)
	# hs <- projectRaster(hs, crs=CRS(madEaProj))

	# # occurrences
	# varecia <- read.csv('./Data/Varecia/03 Added Ratsimbazafy and Parsed Year.csv')
	# varecia <- varecia[varecia$presAbs == 1, ]
	# varecia <- SpatialPointsDataFrame(varecia[ , longLat], data=varecia, proj4=getCRS('wgs84', TRUE))
	# varecia <- sp::spTransform(varecia, CRS(madEaProj))
	
	# # ancillary geo data
	# load('./Study Region & Masks/UTM 38S 30-m Resolution/Eastern Humid Forest Polygon Buffer.RData')
	# load('./Study Region & Masks/UTM 38S 30-m Resolution/Eastern Humid Forest Polygon.RData')
	# load('./Study Region & Masks/UTM 38S 30-m Resolution/Madagascar from GADM 3.6.RData')
	# load('./Data/Protected Areas/WDPA_Sept2018_MDG-shapefile-polygons-onlyTerrestrial.RData')
	
	# # pas <- crop(pas, humidForestBuffer_utm38s)

	# pch <- ifelse(varecia$species == 'Varecia variegata', 21, 24)
	# bg <- ifelse(varecia$species == 'Varecia variegata', 'white', 'red')

	# grays <- paste0('gray', 0:100)
	
	# png(paste0(outDir, '/Occurrence Map for Varecia - Both Species.png'), width=600, height=1200, res=300)
		
		# par(mar=0.5 * c(1, 1, 1, 1), oma=c(1, 1, 1, 1), lwd=0.6)
		
		# plot(humidForest_utm38s, ann=FALSE)
		
		# for (i in seq_along(grays)) grays[i] <- alpha(grays[i], 0.4)
		# plot(madagascar_utm38s)
		# plot(hs, add=TRUE, col=grays, legend=FALSE)
		# plot(madagascar_utm38s, add=TRUE)
		# plot(humidForest_utm38s, add=TRUE, border=alpha('chartreuse', 0.4), col=alpha('chartreuse', 0.1), lwd=1.2)
		# plot(humidForestBuffer_utm38s, add=TRUE, border='darkred', col=NA, lwd=1.2)
		# plot(pas, add=TRUE, col=alpha('blue', 0.20), border='blue')
		# points(varecia, pch=pch, bg=bg, cex=0.5)
		
		# legend('topleft', inset=c(0, 0.15), bty='n', legend=c('V. variegata', 'V. rubra', 'Humid forest', 'Study region', 'Protected'), pch=c(21, 24, NA, NA, NA), col=c('black', 'black', NA, NA, NA), pt.bg=c('white', 'red', NA, NA, NA), border=c(NA, NA, 'chartreuse', 'darkred', 'blue'), fill=c(NA, NA, 'darkseagreen1', NA, alpha('blue', 0.2)), cex=0.5)

		# title(sub=date(), cex.sub=0.3, line=-0, xpd=NA)
		
	# dev.off()
	
# say('#############################################')
# say('### create display of Madagascar with PAs ###')
# say('#############################################')	

	# # no occurrence records in this map!

	# outDir <- './Figures & Tables/'
	# dirCreate(outDir)

	# # hillshading
	# elev <- raster('./Data/Topography - WORLDCLIM Ver 1pt4 Rel 3/elevation_wgs84.tif')
	# slope <- terrain(elev, 'slope')
	# aspect <- terrain(elev, 'aspect')
	
	# hs <- hillShade(slope, aspect, angle=315)
	# hs <- projectRaster(hs, crs=CRS(madEaProj))

	# # ancillary geo data
	# load('./Study Region & Masks/UTM 38S 30-m Resolution/Madagascar from GADM 3.6.RData')
	# load('./Data/Protected Areas/WDPA_Sept2018_MDG-shapefile-polygons-onlyTerrestrial.RData')
	
	# grays <- paste0('gray', 0:100)
	
	# png(paste0(outDir, '/Madagascar with PAs.png'), width=800, height=1200, res=300)
		
		# par(mar=0.1 * c(1, 1, 1, 1), oma=0.1 * c(1, 1, 1, 1))
		
		# for (i in seq_along(grays)) grays[i] <- alpha(grays[i], 0.4)
		# plot(madagascar_utm38s, border=NA, ann=FALSE)
		# plot(hs, col=grays, legend=FALSE)
		# plot(pas, add=TRUE, col=alpha('blue', 0.45), border=NA)
		# title(sub=date(), cex.sub=0.3, line=0, xpd=NA)
		
	# dev.off()
	
# say('########################################################')
# say('### match occurrence records with environmental data ###')
# say('########################################################')

	# varecia <- read.csv('./Data/Varecia/03 Added Ratsimbazafy and Parsed Year.csv')

	# # predictors in WGS84 30 arcsec resolution
	# preds_wgs84 <- stackMadCurrent_wgs84()
	# varecia_wgs84 <- SpatialPoints(varecia[ , longLat], getCRS('wgs84', TRUE))
	# env <- raster::extract(preds_wgs84, varecia_wgs84)
	# env <- as.data.frame(env)

	# varecia <- cbind(varecia, env)
	
	# # predictors in UTM38S 30-m resolution
	# preds_utm38s <- stackMadCurrent_utm38s()
	# varecia_utm38s <- sp::spTransform(varecia_wgs84, CRS(madEaProj))
	# env <- raster::extract(preds_utm38s, varecia_utm38s)
	# env <- as.data.frame(env)
	
	# varecia <- cbind(varecia, env)

	# coords <- coordinates(varecia_utm38s)
	# coords <- as.data.frame(coords)
	# names(coords) <- c('longUtm38s', 'latUtm38s')
	
	# varecia <- insertCol(coords, into=varecia, at='longWgs84', before=FALSE)
	
	# varecia$forest2000_utm38s[is.na(varecia$forest2000_utm38s)] <- 0
	# varecia$forest2005_utm38s[is.na(varecia$forest2005_utm38s)] <- 0
	# varecia$forest2010_utm38s[is.na(varecia$forest2010_utm38s)] <- 0
	# varecia$forest2014_utm38s[is.na(varecia$forest2014_utm38s)] <- 0
	
	# varecia$forestFragClass2000_utm38s <- as.factor(varecia$forestFragClass2000_utm38s)
	# varecia$forestFragClass2005_utm38s <- as.factor(varecia$forestFragClass2005_utm38s)
	# varecia$forestFragClass2010_utm38s <- as.factor(varecia$forestFragClass2010_utm38s)
	# varecia$forestFragClass2014_utm38s <- as.factor(varecia$forestFragClass2014_utm38s)
	
	# ### fil NAs: for sites in water use values from nearest terrestrial cell
	# ########################################################################
	
		# this <- which(varecia$species == 'Varecia rubra' & round(varecia$latWgs84, 5) == -15.67847 & varecia$longWgs84 == 49.9579 & varecia$provider == 'Borgerson' & varecia$origDate == '2011-2012')

			# new <- orig <- varecia[this, ]
			# new$longWgs84 <- 49.958594861
			# new <- raster::extract(preds_wgs84, new[ , longLat])
			# new <- as.data.frame(new)
			
			# for (name in names(new)) if (is.na(orig[ , name])) orig[ , name] <- new[ , name]
			# varecia[this, ] <- orig
		
		# this <- which(varecia$species == 'Varecia rubra' & round(varecia$latWgs84, 5) == -15.67392 & round(varecia$longWgs84, 5) == 49.95717 & varecia$provider == 'Borgerson' & varecia$origDate == '2011-2012')
		
			# new <- orig <- varecia[this, ]
			# new$longWgs84 <- 49.958594861
			# new <- raster::extract(preds_wgs84, new[ , longLat])
			# new <- as.data.frame(new)
			
			# for (name in names(new)) if (is.na(orig[ , name])) orig[ , name] <- new[ , name]
			# varecia[this, ] <- orig

		# this <- which(varecia$species == 'Varecia variegata' & round(varecia$latWgs84, 5) == -15.49833 & round(varecia$longWgs84, 5) == 49.76567 & varecia$provider == 'EEL' & varecia$origDate == '10/12/2004')
		
			# new <- orig <- varecia[this, ]
			# new$longWgs84 <- 49.766970898
			# new <- raster::extract(preds_wgs84, new[ , longLat])
			# new <- as.data.frame(new)
			
			# for (name in names(new)) if (is.na(orig[ , name])) orig[ , name] <- new[ , name]
			# varecia[this, ] <- orig
		
		# this <- which(varecia$species == 'Varecia variegata' & round(varecia$latWgs84, 5) == -15.49308 & round(varecia$longWgs84, 5) == 49.76058 & varecia$provider == 'EEL' & varecia$origDate == '10/16/2004')
		
			# new <- orig <- varecia[this, ]
			# new$latWgs84 <- -15.491430824
			# new <- raster::extract(preds_wgs84, new[ , longLat])
			# new <- as.data.frame(new)
			
			# for (name in names(new)) if (is.na(orig[ , name])) orig[ , name] <- new[ , name]
			# varecia[this, ] <- orig

	# write.csv(varecia, './Data/Varecia/04 Varecia Occurrences Matched with Environmental Data.csv', row.names=FALSE)

# say('###############################')
# say('### select background sites ###')
# say('###############################')

	# elev <- raster('./Data/Topography - WORLDCLIM Ver 1pt4 Rel 3/elevation_wgs84.tif')
	
	# sites <- randomPoints(elev, 11000)
	# sites <- as.data.frame(sites)
	# names(sites) <- longLat
	
	# sites_wgs84 <- SpatialPoints(sites, getCRS('wgs84', TRUE))

	# # predictors in WGS84 30 arcsec resolution
	# preds_wgs84 <- stackMadCurrent_wgs84()
	
	# env <- raster::extract(preds_wgs84, sites_wgs84)
	# env <- as.data.frame(env)
	
	# env$forest2000_wgs84[is.na(env$forest2000_wgs84)] <- 0
	# env$forest2005_wgs84[is.na(env$forest2005_wgs84)] <- 0
	# env$forest2010_wgs84[is.na(env$forest2010_wgs84)] <- 0
	# env$forest2014_wgs84[is.na(env$forest2014_wgs84)] <- 0
	
	# bg <- cbind(sites, env)
	
	# # predictors in UTM38S 30-m resolution
	# preds_utm38s <- stackMadCurrent_utm38s()
	# sites_utm38s <- sp::spTransform(sites_wgs84, CRS(madEaProj))
	# env <- raster::extract(preds_utm38s, sites_utm38s)
	# env <- as.data.frame(env)
	
	# bg <- cbind(bg, env)

	# coords <- coordinates(sites_utm38s)
	# coords <- as.data.frame(coords)
	# names(coords) <- c('longUtm38s', 'latUtm38s')
	
	# bg <- insertCol(coords, into=bg, at='longWgs84', before=FALSE)
	
	# bg$forest2000_utm38s[is.na(bg$forest2000_utm38s)] <- 0
	# bg$forest2005_utm38s[is.na(bg$forest2005_utm38s)] <- 0
	# bg$forest2010_utm38s[is.na(bg$forest2010_utm38s)] <- 0
	# bg$forest2014_utm38s[is.na(bg$forest2014_utm38s)] <- 0
	
	# bg$forestFragClass2000_utm38s[is.na(bg$forestFragClass2000_utm38s)] <- 0
	# bg$forestFragClass2005_utm38s[is.na(bg$forestFragClass2005_utm38s)] <- 0
	# bg$forestFragClass2010_utm38s[is.na(bg$forestFragClass2010_utm38s)] <- 0
	# bg$forestFragClass2014_utm38s[is.na(bg$forestFragClass2014_utm38s)] <- 0

	# bg$forestFragClass2000_utm38s <- as.factor(bg$forestFragClass2000_utm38s)
	# bg$forestFragClass2005_utm38s <- as.factor(bg$forestFragClass2005_utm38s)
	# bg$forestFragClass2010_utm38s <- as.factor(bg$forestFragClass2010_utm38s)
	# bg$forestFragClass2014_utm38s <- as.factor(bg$forestFragClass2014_utm38s)

	# bg$forestFragConnect2000_utm38s[is.na(bg$forestFragConnect2000_utm38s)] <- 0
	# bg$forestFragConnect2005_utm38s[is.na(bg$forestFragConnect2005_utm38s)] <- 0
	# bg$forestFragConnect2010_utm38s[is.na(bg$forestFragConnect2010_utm38s)] <- 0
	# bg$forestFragConnect2014_utm38s[is.na(bg$forestFragConnect2014_utm38s)] <- 0

	# bg$forestFragDensity2000_utm38s[is.na(bg$forestFragDensity2000_utm38s)] <- 0
	# bg$forestFragDensity2005_utm38s[is.na(bg$forestFragDensity2005_utm38s)] <- 0
	# bg$forestFragDensity2010_utm38s[is.na(bg$forestFragDensity2010_utm38s)] <- 0
	# bg$forestFragDensity2014_utm38s[is.na(bg$forestFragDensity2014_utm38s)] <- 0

	# # remove rows with any NAs
	# nas <- naRows(bg)
	# if (length(bg) > 0) bg <- bg[-nas, ]
	
	# bg <- bg[1:10000, ]
	# rownames(bg) <- 1:10000

	# dirCreate('./Ecological Niche Models')
	# save(bg, file='./Ecological Niche Models/Random Background Sites from across Madagascar with Current Environmental Data.RData')
	
# say('##############################################################')
# say('### calculate spatial autocorrelation between survey sites ###')
# say('##############################################################')
	
	# say('I want to use inverse-p weighting applied to survey sites to correct for spatial sampling bias. I will assume that two survey sites at the exact same location each have a weight of 1/2, three each have a weight of 1/3, and so on. At the other end of the spectrum survey sites that are far enough away should each have a weight of 1. I will define "far enough away" as the distance at which the proportion of pairwise observed distances falls below the upper 95th quantile of the distribution of pairwise distances from randomly located sites (one-tail). I will draw a number of randomly located sites in each iteration so it is equal to the total number of survey sites.', breaks=80)

	# varecia <- read.csv('./Data/Varecia/04 Varecia Occurrences Matched with Environmental Data.csv')
	# load('./Ecological Niche Models/Random Background Sites from across Madagascar with Current Environmental Data.RData')

	# # remove duplicates (often survey sites for one species are also counted as surveys for the other)
	# dups <- integer()
	# for (i in 1:nrow(varecia)) {
	
		# this <- varecia[i, ]
		# matches <- which(this$longWgs84 %==% varecia$longWgs84 &
			# this$latWgs84 %==% varecia$latWgs84 &
			# this$provider == varecia$provider &
			# this$origDate == varecia$origDate
		# )
		
		# matches <- matches[-which(matches <= i)]
		# dups <- c(dups, matches)
	
	# }
	
	# vareciaNoDups <- varecia[-dups, ]
	
	# # observed inter-point distances
	# vareciaNoDups <- SpatialPoints(vareciaNoDups[ , longLat], getCRS('wgs84', TRUE))
	# vareciaDists <- distm(vareciaNoDups)
	# vareciaDists <- c(vareciaDists)
	# vareciaDists <- vareciaDists / 1000
	
	# distStep <- 20
	# maxUpperDist <- distStep * ceiling(max(vareciaDists, na.rm=TRUE) / distStep)

	# breaks <- matrix(
		# c(
			# seq(0, maxUpperDist, by=0.5 * distStep),
			# seq(0, maxUpperDist, by=0.5 * distStep) + distStep
		# ),
		# ncol=2,
		# byrow=FALSE
	# )

	# distDistrib <- histOverlap(vareciaDists, breaks=breaks, right=FALSE, graph=FALSE)
	
	# # generate random inter-point distances using same number of points as observed
	
	# randDists <- pointDist(bg[ , longLat])
	# randDists <- randDists / 1000
	# randIndex <- 1:nrow(bg)
	
	# # iterate
	# for (i in 1:100) {
		
		# say(i)
	
		# theseRand <- randDists[sample(randIndex, length(vareciaNoDups), replace=TRUE), ]
		# theseRand <- c(theseRand)
		# randDistDistrib <- histOverlap(theseRand, breaks=breaks, right=FALSE, graph=FALSE)
		# distDistrib <- cbind(distDistrib, randDistDistrib[ , 'proportion'])
		# colnames(distDistrib)[ncol(distDistrib)] <- paste0('rand', prefix(i, 3))
	
	# }
	
	# # plot
	# dirCreate('./Figures & Tables/Spatial Autocorrelation between Survey Sites')
	# png('./Figures & Tables/Spatial Autocorrelation between Survey Sites/Comparison of Distances between Survey Sites and Random Sites.png', width=1200, height=600, res=150)
		
		# ylim <- c(0, max(c(distDistrib[ , 'proportion'], distDistrib[ , grepl(colnames(distDistrib), pattern='rand')])))
		# ats <- rowMeans(distDistrib[ , 1:2])
		# plot(ats, ats, ylim=ylim, col='white', main='Distances between Survey Sites and Random Sites', xlab='Distance (km)', ylab='Proportion of Sites')
		
		# sigSoFar <- TRUE # flips on first time observed proportion of distances falls in 95% CI of randomized distribution of ranges
		# countArrows <- 0
		# for (i in 1:nrow(distDistrib)) {
		
			# at <- ats[i]
			
			# randProp <- c(distDistrib[i , grepl(colnames(distDistrib), pattern='rand')])
			# randPropQuants <- quantile(randProp, c(0.05, 0.5, 0.95))
			# points(at, randPropQuants[3], pch=5, col='blue')
			
			# obsProp <- distDistrib[i , 'proportion']
			# insig <- (obsProp <= randPropQuants[3])
			# pch <- if (insig) { 1 } else { 16 }
			# points(at, obsProp, pch=pch, cex=1)
			
			# if (insig & sigSoFar) {
				# countArrows <- countArrows + 1
				# if (countArrows == 1) {
					# charDist <- at
					# arrows(x0=at, x1=at, y0=-0.19 * ylim[2], y1=-0.0025, col='red', angle=15, length=0.1, xpd=NA, lwd=1.5)
				# }
			# }
			
		# }
		
		# legend('topright', inset=0.01, bty='n', legend=c('Observed (not significant)', 'Observed (significant)', 'Random (Upper 95% CI)'), pch=c(1, 16, 5), col=c('black', 'black', 'blue'))
		
	# dev.off()
		
	# out <- data.frame(characteristicDistanceOfSurveys_km = charDist)
	# write.csv(out, './Figures & Tables/Spatial Autocorrelation between Survey Sites/Characteristic Distance at which Presence of Survey Sites are No Longer Correlated.csv', row.names=FALSE)

# say('##################################################################')
# say('### calculate survey site weights to correct for sampling bias ###')
# say('##################################################################')

	# varecia <- read.csv('./Data/Varecia/04 Varecia Occurrences Matched with Environmental Data.csv')
	# charDist <- read.csv('./Figures & Tables/Spatial Autocorrelation between Survey Sites/Characteristic Distance at which Presence of Survey Sites are No Longer Correlated.csv')
	# charDist <- charDist$characteristicDistanceOfSurveys_km
	
	# # calculate pairwise distances
	# vareciaDists <- pointDist(varecia[ , longLat])
	# diag(vareciaDists) <- NA

	# # calculate weights
	# weight <- rep(NA, nrow(varecia))
	# for (i in seq_along(weight)) {

		# # get only sites within characteristic distance that are for this species
		# neighIndex <- which(vareciaDists[i , ] < charDist & varecia$species[i] == varecia$species)
		# neighDists <- vareciaDists[i, neighIndex]

		# if (length(neighDists) == 0) {
			# weight[i] <- 1
		# } else {

			# numNeighs <- length(neighDists)
			# effectiveNeighs <- sum(1 - (charDist - neighDists) / charDist)
			# weight[i] <- (1 + effectiveNeighs) / (1 + numNeighs)
			
		# }
		
	# }
	
	# varecia$weight <- weight
	# write.csv(varecia, './Data/Varecia/05 Varecia Occurrences with Weights Based on Neighbor Proximity.csv', row.names=FALSE)
		
# say('#######################################################')
# say('### collate occurrence records and background sites ###')
# say('#######################################################')

	# load('./Ecological Niche Models/Random Background Sites from across Madagascar with Current Environmental Data.RData')
	# varecia <- read.csv('./Data/Varecia/05 Varecia Occurrences with Weights Based on Neighbor Proximity.csv')
	
	# bg$weight <- 1
	
	# preds_wgs84 <- stackMadCurrent_wgs84()
	# preds_utm38s <- stackMadCurrent_utm38s()
	
	# for (taxon in taxa) {
	
		# say(taxon)
	
		# ### assign geographic folds to presences
	
		# pres <- if (taxon == 'genus') {
			# varecia[varecia$presAbs == 1, ]
		# } else {
			# varecia[varecia$species == paste('Varecia', taxon) & varecia$presAbs == 1, ]
		# }

		# if (taxon == 'genus') {
		
			# presFold <- rep(NA, nrow(pres))
		
			# presFold[pres$latWgs84 > -16.06] <- 1
			# presFold[pres$latWgs84 <= -16.06 & pres$latWgs84 > -19.61] <- 2
			# presFold[pres$latWgs84 <= -19.61] <- 3
		
		# } else if (taxon == 'variegata') {
		
			# presFold <- rep(NA, nrow(pres))
		
			# presFold[pres$latWgs84 > -17.20] <- 1
			# presFold[pres$latWgs84 <= -17.20 & pres$latWgs84 > -19.67] <- 2
			# presFold[pres$latWgs84 <= -19.67] <- 3

		# } else if (taxon == 'rubra') {
		
			# presFold <- geoFold(pres[ , longLat], k=3, minIn=15, swaps=0)
		
		# }
		
		# ### assign background sites to fold of nearest presence site
		
		# dists <- pointDist(pres[ , longLat], bg[ , longLat])
		# closestSiteToBg <- apply(dists, 2, which.min)
		
		# bgFold <- presFold[closestSiteToBg]
		
		# foldName <- paste0('fold', capIt(taxon))
		
		# pres$DUMMY <- presFold
		# names(pres)[ncol(pres)] <- foldName
		
		# bg$DUMMY <- bgFold
		# names(bg)[ncol(bg)] <- foldName
		
		# bg$year <- NA
		
		# # combine presences and background sites
		# presBg <- c(rep(1, nrow(pres)), rep(0, nrow(bg)))
		
		# taxonData <- rbind(
			# pres[ , c('longWgs84', 'latWgs84', 'year', names(preds_wgs84), names(preds_utm38s), 'weight', foldName)],
			# bg[ , c('longWgs84', 'latWgs84', 'year', names(preds_wgs84), names(preds_utm38s), 'weight', foldName)]
		# )
		
		# taxonData <- cbind(presBg, taxonData)
		
		# ### define "forestCover" and "forestFragClass" variables as forest cover or fragmentation class of year closest to the year of survey
		# # use 2005 for background sites and survey sites with no date
		
		# taxonData$forestCover_wgs84 <- NA
		# taxonData$forestCover_utm38s <- NA
		# taxonData$forestFragClass_utm38s <- NA
		# forestYears <- c(2000, 2005, 2010, 2014)

		# # assign forest of survey sites WITH census year
		# yearNotNa <- which(!is.na(taxonData$year))
		
		# if (length(yearNotNa) > 0) {
			
			# for (thisRecord in yearNotNa) {
			
				# yearDiff <- which.min(abs(forestYears - taxonData$year[thisRecord]))
				# nearestYear <- forestYears[yearDiff]
				
				# taxonData$forestCover_wgs84[thisRecord] <- taxonData[thisRecord, paste0('forest', nearestYear, '_wgs84')]
				# taxonData$forestCover_utm38s[thisRecord] <- taxonData[thisRecord, paste0('forest', nearestYear, '_utm38s')]
			
				# taxonData$forestFragClass_utm38s[thisRecord] <- taxonData[thisRecord, paste0('forestFragClass', nearestYear, '_utm38s')]
			
			# }
			
		# }

		# # assign forest of survey sites WITHOUT census year
		# medianYear <- median(taxonData$year, na.rm=TRUE)
		# yearDiff <- which.min(abs(forestYears - medianYear))
		# nearestYear <- forestYears[yearDiff]
		
		# yearNa <- which(is.na(taxonData$year))
		# taxonData$forestCover_wgs84[yearNa] <- taxonData[yearNa , paste0('forest', nearestYear, '_wgs84')]
		# taxonData$forestCover_utm38s[yearNa] <- taxonData[yearNa , paste0('forest', nearestYear, '_utm38s')]
		
		# taxonData$forestFragClass_utm38s[yearNa] <- taxonData[yearNa , paste0('forestFragClass', nearestYear, '_utm38s')]
		
		# save(taxonData, file=paste0('./Ecological Niche Models/Collated Presence and Background Point Data for Varecia ', toupper(taxon), '.RData'))
		
	# } # next taxon

# say('#####################################')
# say('### train ecological niche models ###')
# say('#####################################')	
	
	# for (taxon in taxa) {
	
		# say(taxon, level=1)
		
		# load(paste0('./Ecological Niche Models/Collated Presence and Background Point Data for Varecia ', toupper(taxon), '.RData'))
		
		# taxonData$forestFragClass_utm38s <- as.factor(taxonData$forestFragClass_utm38s)
		# taxonData[ , paste0('bio', prefix(bioclims, 2))] <- apply(taxonData[ , paste0('bio', prefix(bioclims, 2))], 2, stretchMinMax)
		
		# foldName <- paste0('fold', capIt(taxon))
		
		# ### by geo-fold
		# ###############
		
		# dirCreate('./Ecological Niche Models')
		
		# for (k in 1:3) {
		
			# say('k = ', k, ' =========================================================')
		
			# trainData <- taxonData[taxonData[ , foldName] != k, ]
			# trainWeight <- balanceWeights(trainData)

			# model <- trainGlm(resp='presBg', preds=predictors, data=trainData, w=trainWeight, verboten=c('I(forestFragClass_utm38s^2)'), method='brglmFit', verbose=TRUE)
		
			# save(model, file=paste0('./Ecological Niche Models/GLM Niche Model for Varecia ', taxon, ' for fold ', k, '.RData'))
		
			# say('', post=2)
		
			# ### NS
			# model <- trainNs(data=trainData, resp='presBg', preds=predictors, w=trainWeight, verbose=TRUE)
			
			# save(model, file=paste0('./Ecological Niche Models/NS Niche Model for Varecia ', taxon, ' for fold ', k, '.RData'))
		
			# say('', post=4)
		
		# }

		# ### all-sites models
		# ####################
		
		# trainData <- taxonData
		# trainWeight <- balanceWeights(trainData)

		# ### GLM
		# model <- trainGlm(data=trainData, resp='presBg', preds=predictors, w=trainWeight, verboten=c('I(forestFragClass_utm38s^2)'), method='brglmFit', verbose=TRUE)
		
		# save(model, file=paste0('./Ecological Niche Models/GLM Niche Model for Varecia ', taxon, ' all-sites model.RData'))
	
		# ### NS
		# model <- trainNs(data=trainData, resp='presBg', preds=predictors, w=trainWeight, verbose=TRUE)
		
		# save(model, file=paste0('./Ecological Niche Models/NS Niche Model for Varecia ', taxon, ' all-sites model.RData'))
	
	# } # next taxon
	
# say('####################################################################################')
# say('### select and collate future climate data for projecting ecological niche model ###')
# say('####################################################################################')

	# mask_wgs84 <- raster('./Study Region & Masks/WGS84 30-arcsec Resolution/mask_wgs84.tif')

	# ### calculate differences between current and future climate
	# ############################################################
	
	# say('Selecting ESMs based on which ones give the most extreme values of the selected BIOCLIMs using the 2061-2080 period under RCP 8.5.', breaks=80)
	
	# say('Calculating differences between current and future climate for select variables...', pre=2)
	
	# gcmInfo <- read.csv('F:/ecology/Climate/WORLDCLIM Ver 1pt4 Rel 3/GCM Names, Scenarios Modeled, and File Codes.csv')
	# gcmsUsed <- c(1, which(gcmInfo$useToEnsemble))

	# diffsByGcm <- data.frame()
	
	# load('./Ecological Niche Models/Collated Presence and Background Point Data for Varecia GENUS.RData')
	# genus <- taxonData[taxonData$presBg == 1, longLat]
	# load('./Ecological Niche Models/Collated Presence and Background Point Data for Varecia VARIEGATA.RData')
	# variegata <- taxonData[taxonData$presBg == 1, longLat]
	# load('./Ecological Niche Models/Collated Presence and Background Point Data for Varecia RUBRA.RData')
	# rubra <- taxonData[taxonData$presBg == 1, longLat]
	
	# variegata <- elimCellDups(variegata, mask_wgs84)
	# rubra <- elimCellDups(rubra, mask_wgs84)
	# genus <- elimCellDups(genus, mask_wgs84)
	
	# for (thisBio in bioclims) {
		
		# say(thisBio, pre=1, post=0)
		
		# currentFile <- paste0('./Data/Climate - WORLDCLIM Ver 2 Rel 1 Jun 2016 (1970-2000)/bio', prefix(thisBio, 2), '.tif')
		# currentClimate <- raster(currentFile)
		
		# for (countGcm in gcmsUsed) {
		
			# gcm <- gcmInfo$gcm[countGcm]
			# gcmCode <- tolower(gcmInfo$code[countGcm])

			# say(gcm, post=0)
			
			# futureFile <- paste0('F:/ecology/Climate/WORLDCLIM Ver 1pt4 Rel 3/30 arcsec/RCP8pt5 - 30 arcsec/2061 to 2080 - 30 arcsec - RCP8pt5 - ', gcm, '/World/', gcmCode, '85bi70', prefix(thisBio, 2), '.tif')
			# futureClimate <- raster(futureFile)
			# futureClimate <- crop(futureClimate, mask_wgs84)
			
			# if (thisBio %in% c(4, 10)) futureClimate <- futureClimate / 10
			
			# diff <- futureClimate - currentClimate
			# deltaAllOfMad <- cellStats(diff, 'mean')
			
			# deltaVariegata <- mean(raster::extract(diff, variegata), na.rm=TRUE)
			# deltaRubra <- mean(raster::extract(diff, rubra), na.rm=TRUE)
			# deltaGenus <- mean(raster::extract(diff, genus), na.rm=TRUE)
			
			# diffsByGcm <- rbind(
				# diffsByGcm,
				# data.frame(
					# gcm=gcm,
					# bio=paste0('bio', prefix(thisBio, 2)),
					# deltaAllOfMad=deltaAllOfMad,
					# deltaVariegata=deltaVariegata,
					# deltaRubra=deltaRubra,
					# deltaGenus=deltaGenus
				# )
			# )
	
		# } # next BIOCLIM
		
	# } # next ESM
	
	# dirCreate('./Figures & Tables/Climate Model Selection')
	# write.csv(diffsByGcm, './Figures & Tables/Climate Model Selection/Differences between Current and Future Values for each BIOCLIM by GCM.csv', row.names=FALSE)

	# ### select ESMs with greatest differences between themselves across occupied cells of genus
	# ###########################################################################################
	
	# say('Identifying ESM with most extreme values from present...', pre=2)
	
	# useTheseGcms <- character()
	
	# for (bio in bioclims) {
	
		# bioclim <- paste0('bio', prefix(bio, 2))
		# candidates <- diffsByGcm[diffsByGcm$bio == bioclim, ]
		# if (any(candidates$deltaGenus < 0)) useTheseGcms <- c(useTheseGcms, candidates$gcm[which.min(candidates$deltaGenus)])
		# if (any(candidates$deltaGenus > 0)) useTheseGcms <- c(useTheseGcms, candidates$gcm[which.max(candidates$deltaGenus)])
		
	# }
		
	# useTheseGcms <- sort(unique(useTheseGcms))
				
	# say('GCMs with most extreme values: ', paste(useTheseGcms, collapse=' '), pre=1)

	# out <- data.frame(gcm=useTheseGcms, col=c('red', 'orange', 'green', 'blue', 'purple'))
	
	# write.csv(out, './Figures & Tables/Climate Model Selection/ESMs with the Greatest Difference between Current and Future Values for each BIOCLIM.csv', row.names=FALSE)
	
	# ###################################
	# ### collate future climate data ###
	# ###################################
	
	# say('Collating future climate data...', pre=2)
	
	# for (gcm in useTheseGcms) {
	
		# gcmCode <- tolower(gcmInfo$code[gcmInfo$gcm == gcm])
	
		# for (period in periods) {
			
			# periodCode <- if (period == '2041 to 2060') {
				# 50
			# } else {
				# 70
			# }
			
			# for (rcp in rcps) {
	
				# say(gcm, ' ', period, ' ', rcp, post=0)

				# rcpCode <- if (rcp == '2pt6') {
					# '26'
				# } else if (rcp == '4pt5') {
					# '45'
				# } else if (rcp == '6pt0') {
					# '60'
				# } else if (rcp == '8pt5') {
					# '85'
				# }
				
				# outDir <- paste0('./Data/Climate - WORLDCLIM Ver 2 Rel 1 Jun 2016 (', period, ')/', gcm, ' RCP', rcp)
				# dirCreate(outDir)
	
				# for (thisBio in c(bioclims)) {
				
					# say(thisBio, post=0)
				
					# rastFile <- paste0('F:/ecology/Climate/WORLDCLIM Ver 1pt4 Rel 3/30 arcsec/RCP', rcp, ' - 30 arcsec/', period, ' - 30 arcsec - RCP', rcp, ' - ', gcm, '/World/', gcmCode, rcpCode, 'bi', periodCode, prefix(thisBio, 2), '.tif')
				
					# rast <- raster(rastFile)
					
					# rast <- crop(rast, mask_wgs84)
					# names(rast) <- paste0('bio', prefix(thisBio, 2))
					# writeRaster(rast, paste0(outDir, '/bio', prefix(thisBio, 2)))
			
				# }
				
				# say('')
				
			# } # rcp
			
		# } # period
		
	# } # gcm

# say('############################################################')
# say('### visualize response curves of ecological niche models ###')
# say('############################################################')
	
	# fragClasses <- c(0, 1, 2, 3, 4, 6)
	
	# fragCols <- c('gray10', 'blue', 'lightblue', 'yellow', 'orange', 'forestgreen')
	# names(fragCols) <- fragClasses
	
	# mask_wgs84 <- raster('./Study Region & Masks/WGS84 30-arcsec Resolution/mask_wgs84.tif')
	
	# extremeGcms <- read.csv('./Figures & Tables/Climate Model Selection/ESMs with the Greatest Difference between Current and Future Values for each BIOCLIM.csv')
	# gcmDiffs <- read.csv('./Figures & Tables/Climate Model Selection/Differences between Current and Future Values for each BIOCLIM by GCM.csv')
	
	# outDir <- './Figures & Tables/Ecological Niche Models - Responses to Variables'
	# dirCreate(outDir)
	
	# for (taxon in taxa) {
		
		# load(paste0('./Ecological Niche Models/Collated Presence and Background Point Data for Varecia ', toupper(taxon), '.RData'))

		# taxonDataElimCellDups <- elimCellDups(taxonData, mask_wgs84, longLat=longLat)
		
		# rescaledTaxonData <- taxonData
		# rescaledTaxonData[ , paste0('bio', prefix(bioclims, 2))] <- apply(rescaledTaxonData[ , paste0('bio', prefix(bioclims, 2))], 2, stretchMinMax)
		
		# # for (algo in c('glm', 'ns')) {
		# for (algo in c('glm')) {
		
			# say(taxon, ' ', algo)
		
			# load(paste0('./Ecological Niche Models/', toupper(algo), ' Niche Model for Varecia ', taxon, ' all-sites model.RData'))
			
			# png(paste0(outDir, '/Responses of Varecia ', taxon, ' to Predictors Based on ', toupper(algo), ' Niche Model.png'), width=1800, height=1200, res=300)
			
			# par(mfcol=c(2, 3), mar=0.5 * c(5, 5, 5, 1), oma=c(1, 2, 1, 1), cex.main=0.9, cex.axis=0.6, cex.lab=0.7, mgp=c(3, 0.5, 0), tck=-0.05)
			
			# for (predictor in predictors) {
			
				# if (predictor != 'forestFragClass_utm38s') {

					# if (exists('predictToData', inherits=FALSE)) rm(predictToData)

					# # create data frame with all bioclims but focal one at median across presence sites
					# for (thisBio in paste0('bio', prefix(bioclims, 2))) {
						# medianAtPres <- median(rescaledTaxonData[rescaledTaxonData$presBg==1, thisBio])
						# col <- data.frame(rep(medianAtPres, 500))
						# col[1, 1] <- col[1, 1] + 0.000000001
						# predictToData <- if (exists('predictToData')) { cbind(predictToData, col) } else { col }
						# names(predictToData)[ncol(predictToData)] <- thisBio
					# }
					
					# # focal bioclim varies across range of all training sites
					# predictToData[ , predictor] <- seq(min(rescaledTaxonData[ , predictor]), max(rescaledTaxonData[ , predictor]), length.out=nrow(predictToData))

					# # unscale focal predictor for plotting
					# minRaw <- min(taxonData[ , predictor])
					# maxRaw <- max(taxonData[ , predictor])
					# x <- predictToData[ , predictor] * (maxRaw - minRaw) + minRaw
					
					# # assume all locations are forest interior class
					# predictToData$forestFragClass_utm38s <- 6
					# predictToData$forestFragClass_utm38s <- as.factor(predictToData$forestFragClass_utm38s)
					# levels(predictToData$forestFragClass_utm38s) <- list(
						# '0' = '0',
						# '1' = '1',
						# '2' = '2',
						# '3' = '3',
						# '4' = '4',
						# '6' = '6'
					# )
						
					# predictions <- predict(model, predictToData, type='response')

					# plot(x, predictions, xlab='', ylab='', main=predictorNice(predictor), type='l', lwd=2, ylim=c(0, 1.05))
						
					# title(ylab='Suitability', line=1.5)
					# title(xlab=predictorNice(predictor, TRUE), line=1.5)
					
					# ### add indicators for future climate
					# currentClimate <- mean(taxonDataElimCellDups[taxonDataElimCellDups$presBg == 1, predictor])

					# usr <- par('usr')
					# posScale <- 0 # y position scale for future climate indicators

					# for (countGcm in 1:nrow(extremeGcms)) {
						
						# gcm <- extremeGcms$gcm[countGcm]
						
						# delta <- gcmDiffs[which(gcmDiffs == gcm & gcmDiffs$bio == predictor), paste0('delta', capIt(taxon))]
						# futureClimate <- currentClimate + delta
						
						# posScale <- posScale + 0.05
						# yPos <- usr[4] - posScale * (usr[4] - usr[3])

						# col <- extremeGcms$col[countGcm]
						# arrows(x0=currentClimate, x1=futureClimate, y0=yPos, y1=yPos, angle=15, length=0.05, xpd=NA, col=col)
								
					# } # next GCM

					# ### add indicators for known presences
					# x <- taxonDataElimCellDups[taxonDataElimCellDups$presBg == 1, predictor]
					# rug(x)

				# } else {
				
					# if (exists('predictToData', inherits=FALSE)) rm(predictToData)

					# for (thisBio in paste0('bio', prefix(bioclims, 2))) {
						# col <- data.frame(rep(median(rescaledTaxonData[rescaledTaxonData$presBg==1, thisBio]), length(fragClasses)))
						# predictToData <- if (exists('predictToData')) { cbind(predictToData, col) } else { col }
						# names(predictToData)[ncol(predictToData)] <- thisBio
					# }

					# predictToData$forestFragClass_utm38s <- fragClasses
					# predictToData$forestFragClass_utm38s <- as.factor(predictToData$forestFragClass_utm38s)
					# levels(predictToData$forestFragClass_utm38s) <- list(
						# '0' = '0',
						# '1' = '1',
						# '2' = '2',
						# '3' = '3',
						# '4' = '4',
						# '6' = '6'
					# )
				
					# predictions <- predict(model, predictToData, type='response')
					
					# names <- c('none', 'patch', 'trans.', 'perf.', 'edge', 'interior')
					# barplot(predictions, names=names, ylab='', main='Forest Fragmentation Class', xlab='Class', col='gray80')
					# title(ylab='Suitability', line=1.5)

				# }
			
			# } # next predictor

			# # legend (its own plot)
			# plot(0, 0, fg='white', ann=FALSE, col.axis='white', col='white')
			# legend('center', bty='n', title='Earth System Model', legend=extremeGcms$gcm, col=extremeGcms$col, lwd=1, cex=0.8)

			# title(main=paste0('Varecia ', taxon, ' Based on ', toupper(algo), ' Niche Model'), outer=TRUE, line=0, cex.main=0.9)
			
			# dev.off()

		# } # next algorithm
	
	# } # next taxon

# say('########################################')
# say('### evaluate ecological niche models ###')
# say('########################################')

	# eval <- data.frame()

	# for (taxon in taxa) {
	
		# load(paste0('./Ecological Niche Models/Collated Presence and Background Point Data for Varecia ', toupper(taxon), '.RData'))
		# taxonData[ , paste0('bio', prefix(bioclims, 2))] <- apply(taxonData[ , paste0('bio', prefix(bioclims, 2))], 2, stretchMinMax)
		
		# for (algo in c('glm', 'ns')) {
				
			# for (k in 1:3) {
	
				# say(taxon, ' ', algo, ' ', k)
				
				# foldName <- paste0('fold', capIt(taxon))
				# trainData <- taxonData[taxonData[ , foldName] != k, ]
				# testData <- taxonData[taxonData[ , foldName] == k, ]
				
				# testWeight <- balanceWeights(testData)
				
				# presWeight <- testWeight[testData$presBg == 1]
				# bgWeight <- testWeight[testData$presBg == 0]

				# # get predictions at test sites
				# load(paste0('./Ecological Niche Models/', toupper(algo), ' Niche Model for Varecia ', taxon, ' for fold ', k, '.RData'))
				
				# predictions <- predict(model, testData, type='response')
				# presPredictions <- predictions[testData$presBg == 1]
				# bgPredictions <- predictions[testData$presBg == 0]
				
				# # calculate performance statistics
				# cbiSpearman <- contBoyce(presPredictions, bgPredictions, presWeight=presWeight, bgWeight=bgWeight, method='spearman', na.rm=TRUE)
				# cbiPearson <- contBoyce(presPredictions, bgPredictions, presWeight=presWeight, bgWeight=bgWeight, method='pearson', na.rm=TRUE)
			
				# thisEval <- data.frame(
					# taxon = taxon,
					# fold = k,
					# algo = algo,
					# trainPres = sum(trainData$presBg == 1),
					# testPres = sum(testData$presBg == 1),
					# trainBg = sum(trainData$presBg == 0),
					# testBg = sum(testData$presBg == 0),
					# cbiSpearman = cbiSpearman,
					# cbiPearson = cbiPearson
				# )
			
				# eval <- rbind(eval, thisEval)
				
			# } # next algorithm
		
		# } # next geo-fold
		
	# } # next taxon
	
	# dirCreate('./Figures & Tables/Ecological Niche Models - Parameterization')
	# write.csv(eval, './Figures & Tables/Ecological Niche Models - Parameterization/Model Evaluation.csv', row.names=FALSE)

	# evalMean <- aggregate(eval, by=list(eval$taxon, eval$algo), mean)
	# evalMean$taxon <- evalMean$algo <- NULL
	# names(evalMean)[1:2] <- c('taxon', 'algo')

	# evalSd <- aggregate(eval, by=list(eval$taxon, eval$algo), sd)
	# evalSd$taxon <- evalSd$algo <- NULL
	# names(evalSd)[1:2] <- c('taxon', 'algo')

	# write.csv(evalMean, './Figures & Tables/Ecological Niche Models - Parameterization/Model Evaluation - Aggregated Using Mean Value.csv', row.names=FALSE)
	# write.csv(evalSd, './Figures & Tables/Ecological Niche Models - Parameterization/Model Evaluation - Aggregated Using SD Value.csv', row.names=FALSE)

# say('################################################')
# say('### predict to current and future conditions ###')
# say('################################################')

	# load('./Ecological Niche Models/Collated Presence and Background Point Data for Varecia GENUS.Rdata')
	# load('./Ecological Niche Models/GLM Niche Model for Varecia genus all-sites model.RData')
	
	# humidForestBufferMask_utm38s <- raster('./Study Region & Masks/UTM 38S 30-m Resolution/humidForestBufferMask_utm38s.tif')
	
	# ### collate CURRENT climate data
	# #################################
		
		# say('### current climate ###')

		# climateCurrent_wgs84 <- stackMadCurrent_wgs84()
		# climateCurrent_wgs84 <- subset(climateCurrent_wgs84, paste0('bio', prefix(bioclims, 2)))
		
		# ### rescale climate data
		
		# for (thisBio in paste0('bio', prefix(bioclims, 2))) {
		
			# minRaw <- min(taxonData[ , thisBio])
			# maxRaw <- max(taxonData[ , thisBio])
			
			# climateCurrent_wgs84[[thisBio]] <- (climateCurrent_wgs84[[thisBio]] - minRaw) / (maxRaw - minRaw)
			
		# }
		
		# ### project climate to UTM38S @ 30-m resolution

		# beginCluster(3)
			# climateCurrent_utm38s <- projectRaster(climateCurrent_wgs84, humidForestBufferMask_utm38s, args=list(method='ngb'))
		# endCluster()

		# names(climateCurrent_utm38s) <- paste0('bio', prefix(bioclims, 2))

	# ### collate FUTURE climate data
	# ###############################
		
		# say('### future climate ###')
		
		# load('./Ecological Niche Models/Collated Presence and Background Point Data for Varecia GENUS.Rdata')

		# gcmInfo <- read.csv('./Figures & Tables/Climate Model Selection/ESMs with the Greatest Difference between Current and Future Values for each BIOCLIM.csv')
		# gcms <- gcmInfo$gcm

		# # by GCM
		# for (gcm in gcms) {
		
			# for (period in periods) {
			
				# periodShort <- if (period == '2041 to 2060') { 2050 } else if (period == '2061 to 2080') { 2070 }
			
				# for (rcp in rcps) {
		
					# say('Processing ', gcm, ' ', period, ' ', rcp, '...')
		
					# climateFuture_wgs84 <- stackMadFuture_wgs84(gcm, period, rcp)
					
					# climateFuture_wgs84[['bio04']] <- climateFuture_wgs84[['bio04']] / 10
					# climateFuture_wgs84[['bio10']] <- climateFuture_wgs84[['bio10']] / 10
					
					# ### rescale climate data
					# for (thisBio in paste0('bio', prefix(bioclims, 2))) {
					
						# minRaw <- min(taxonData[ , thisBio])
						# maxRaw <- max(taxonData[ , thisBio])
						
						# climateFuture_wgs84[[thisBio]] <- (climateFuture_wgs84[[thisBio]] - minRaw) / (maxRaw - minRaw)
						
					# }
					
					# ### project climate to UTM38S @ 30-m resolution
					
					# beginCluster(4)
						# climateFuture_utm38s <- projectRaster(climateFuture_wgs84, humidForestBufferMask_utm38s, args=list(method='ngb'))
					# endCluster()
					
					# climateFuture_utm38s <- humidForestBufferMask_utm38s * climateFuture_utm38s
					
					# names(climateFuture_utm38s) <- paste0('bio', prefix(bioclims, 2))
					
					# gcmName <- gsub(gcm, pattern='-', replacement='')
					# gcmName <- tolower(gcmName)
					# gcmName <- capIt(gcmName)
					# assign(paste0('climateFuture', gcmName, periodShort, 'Rcp', rcp), climateFuture_utm38s)
					
				# } # rcp
				
			# } # period
			
		# } # gcm
		
	# ### collate CURRENT and FUTURE forest fragmentation class
	# #########################################################

		# say('### current and future forest ###')

		# # current forest fragmentation class
		# forestFragClassCurrent_utm38s <- raster('./Data/Forest - Vieilledent et al 2018/forestFragClass2014_utm38s.tif')

		# # future forest fragmentation class assuming protected areas as-is
		# forestFragClass2050_utm38s <- raster('./Deforestation Models/Forest Fragmentation Class 2015-2080 Assuming Pop-Sensitive Loss Amount/forestFragClass2050_utm38s.tif')
		# forestFragClass2070_utm38s <- raster('./Deforestation Models/Forest Fragmentation Class 2015-2080 Assuming Pop-Sensitive Loss Amount/forestFragClass2070_utm38s.tif')
		
		# # future forest fragmentation class assuming protected areas suffer no new deforestation
		# forestFragClass2050PasMasked_utm38s <- raster('./Deforestation Models/Forest Fragmentation Class 2015-2080 Assuming Pop-Sensitive Loss Amount PAs Have 2014 Cover/forestFragClass2050_utm38s.tif')
		# forestFragClass2070PasMasked_utm38s <- raster('./Deforestation Models/Forest Fragmentation Class 2015-2080 Assuming Pop-Sensitive Loss Amount PAs Have 2014 Cover/forestFragClass2070_utm38s.tif')
		
		# names(forestFragClassCurrent_utm38s) <- 'forestFragClass_utm38s'
		# names(forestFragClass2050_utm38s) <- 'forestFragClass_utm38s'
		# names(forestFragClass2070_utm38s) <- 'forestFragClass_utm38s'
		# names(forestFragClass2050PasMasked_utm38s) <- 'forestFragClass_utm38s'
		# names(forestFragClass2070PasMasked_utm38s) <- 'forestFragClass_utm38s'
		
	# ### predict!!!
	# ##############
	
		# dirCreate('./Ecological Niche Models/Prediction Rasters')

		# ### predict to CURRENT climate and CURRENT forest
		# #################################################

		# say('### predicting to current climate and current forest ###')
		
		# preds <- stack(climateCurrent_utm38s, forestFragClassCurrent_utm38s)
		# prediction <- predictVareciaEnm(model, preds)
		# names(prediction) <- 'glmEnm_vareciaGenus_climateCurrent_forest2014_utm38s'

		# writeRaster(prediction, './Ecological Niche Models/Prediction Rasters/glmEnm_vareciaGenus_climateCurrent_forest2014_utm38s', datatype='INT1U')

		# ### predict to CURRENT climate and 2050 forest with defo ANYWHERE
		# #################################################################

		# say('### predicting to current climate and 2050s forest with deforestation anywhere ###')
		
		# preds <- stack(climateCurrent_utm38s, forestFragClass2050_utm38s)
		# prediction <- predictVareciaEnm(model, preds)
		# names(prediction) <- 'glmEnm_vareciaGenus_climateCurrent_forest2050_defoAnywhere_utm38s'

		# writeRaster(prediction, './Ecological Niche Models/Prediction Rasters/glmEnm_vareciaGenus_climateCurrent_forest2050_defoAnywhere_utm38s', datatype='INT1U')
	
		# ### predict to CURRENT climate and 2070 forest with defo ANYWHERE
		# #################################################################

		# say('### predicting to current climate and 2070s forest with deforestation anywhere ###')
		
		# preds <- stack(climateCurrent_utm38s, forestFragClass2070_utm38s)
		# prediction <- predictVareciaEnm(model, preds)
		# names(prediction) <- 'glmEnm_vareciaGenus_climateCurrent_forest2070_defoAnywhere_utm38s'

		# writeRaster(prediction, './Ecological Niche Models/Prediction Rasters/glmEnm_vareciaGenus_climateCurrent_forest2070_defoAnywhere_utm38s', datatype='INT1U')
	
		# ### predict to CURRENT climate and 2050 forest with defo OUTSIDE PAs
		# ####################################################################

		# say('### predicting to current climate and 2050s forest with deforestation outside protected areas ###')
		
		# preds <- stack(climateCurrent_utm38s, forestFragClass2050PasMasked_utm38s)
		# prediction <- predictVareciaEnm(model, preds)
		# names(prediction) <- 'glmEnm_vareciaGenus_climateCurrent_forest2050_defoOutsidePAs_utm38s'

		# writeRaster(prediction, './Ecological Niche Models/Prediction Rasters/glmEnm_vareciaGenus_climateCurrent_forest2050_defoOutsidePAs_utm38s', datatype='INT1U')
	
		# ### predict to CURRENT climate and 2070 forest with defo OUTSIDE PAs
		# ####################################################################

		# say('### predicting to current climate and 2070s forest with deforestation outside protected areas ###')
		
		# preds <- stack(climateCurrent_utm38s, forestFragClass2070PasMasked_utm38s)
		# prediction <- predictVareciaEnm(model, preds)
		# names(prediction) <- 'glmEnm_vareciaGenus_climateCurrent_forest2070_defoOutsidePAs_utm38s'

		# writeRaster(prediction, './Ecological Niche Models/Prediction Rasters/glmEnm_vareciaGenus_climateCurrent_forest2070_defoOutsidePAs_utm38s', datatype='INT1U')
		
		# ### predict to FUTURE climate and CURRENT forest
		# ################################################
		
		# say('### predicting to future climate and current forest ###')

		# gcmInfo <- read.csv('./Figures & Tables/Climate Model Selection/ESMs with the Greatest Difference between Current and Future Values for each BIOCLIM.csv')
		# gcms <- gcmInfo$gcm
		
		# for (gcm in gcms) {
		# # for (gcm in 'BCC-CSM1-1') {
		
			# gcmName <- gsub(gcm, pattern='-', replacement='')
			# gcmName <- tolower(gcmName)
			# gcmName <- capIt(gcmName)

			# for (period in periods) {
				
				# periodShort <- if (period == '2041 to 2060') { 2050 } else if (period == '2061 to 2080') { 2070 }
				
				# for (rcp in rcps) {

					# say(gcm, ' ', period, ' ', rcp)
				
					# thisClimName <- paste0('climateFuture', gcmName, periodShort, 'Rcp', rcp)
					# futureCurrent_utm38s <- get(thisClimName)

					# preds <- stack(futureCurrent_utm38s, forestFragClassCurrent_utm38s)
					# prediction <- predictVareciaEnm(model, preds, 4)
					# name <- paste0('glmEnm_vareciaGenus_climate', periodShort, gcmName, 'Rcp', rcp, '_forest2014_utm38s')
					# names(prediction) <- name

					# writeRaster(prediction, paste0('./Ecological Niche Models/Prediction Rasters/', name), datatype='INT1U')
		
				# } # next RCP
				
			# } # next period
			
		# } # next gcm
		
		# ### predict to CURRENT climate and FUTURE forest with deforestation ANYWHERE
		# ###########################################################################
		
		# say('### predicting to current climate and future forest with deforestation anywhere ###')

		# for (period in periods) {
			
			# say(period)
			
			# periodShort <- if (period == '2041 to 2060') { 2050 } else if (period == '2061 to 2080') { 2070 }

			# forestFragClassFutureName <- paste0('forestFragClass', periodShort, '_utm38s')
			# forestFragClassFuture <- get(forestFragClassFutureName)
			
			# preds <- stack(climateCurrent_utm38s, forestFragClassFuture)
			# prediction <- predictVareciaEnm(model, preds)
			# name <- paste0('glmEnm_vareciaGenus_climateCurrent_forest', periodShort, '_defoAnywhere_utm38s')
			# names(prediction) <- name

			# writeRaster(prediction, paste0('./Ecological Niche Models/Prediction Rasters/', name), datatype='INT1U')
			
		# } # next period

		# ### predict to CURRENT climate and FUTURE forest with deforestation OUTSIDE PROTECTED AREAS
		# ###########################################################################################
		
		# say('### predicting to current climate and future forest with deforestation outside protected areas ###')

		# for (period in periods) {
			
			# say(period)
			
			# periodShort <- if (period == '2041 to 2060') { 2050 } else if (period == '2061 to 2080') { 2070 }

			# forestFragClassFutureName <- paste0('forestFragClass', periodShort, 'PasMasked_utm38s')
			# forestFragClassFuture <- get(forestFragClassFutureName)
			
			# preds <- stack(climateCurrent_utm38s, forestFragClassFuture)
			# prediction <- predictVareciaEnm(model, preds)
			# name <- paste0('glmEnm_vareciaGenus_climateCurrent_forest', periodShort, '_defoOutsidePAs_utm38s')
			# names(prediction) <- name

			# writeRaster(prediction, paste0('./Ecological Niche Models/Prediction Rasters/', name), datatype='INT1U')
			
		# } # next period

		# ### predict to FUTURE climate and FUTURE forest with deforestation ANYWHERE
		# ###########################################################################
		
		# say('### predicting to future climate and future forest with deforestation anywhere ###')
		
		# gcmInfo <- read.csv('./Figures & Tables/Climate Model Selection/ESMs with the Greatest Difference between Current and Future Values for each BIOCLIM.csv')
		# gcms <- gcmInfo$gcm
		
		# for (gcm in gcms) {
		# # for (gcm in 'BCC-CSM1-1') {
		
			# gcmName <- gsub(gcm, pattern='-', replacement='')
			# gcmName <- tolower(gcmName)
			# gcmName <- capIt(gcmName)

			# for (period in periods) {
				
				# periodShort <- if (period == '2041 to 2060') { 2050 } else if (period == '2061 to 2080') { 2070 }

				# forestFragClassFutureName <- paste0('forestFragClass', periodShort, '_utm38s')
				# forestFragClassFuture <- get(forestFragClassFutureName)
				
				# for (rcp in rcps) {

					# say(gcm, ' ', period, ' ', rcp)
				
					# thisClimName <- paste0('climateFuture', gcmName, periodShort, 'Rcp', rcp)
					# futureCurrent_utm38s <- get(thisClimName)

					# preds <- stack(futureCurrent_utm38s, forestFragClassFuture)
					# prediction <- predictVareciaEnm(model, preds, 4)
					# name <- paste0('glmEnm_vareciaGenus_climate', periodShort, gcmName, 'Rcp', rcp, '_forest', periodShort, '_defoAnywhere_utm38s')
					# names(prediction) <- name

					# writeRaster(prediction, paste0('./Ecological Niche Models/Prediction Rasters/', name), datatype='INT1U')
		
				# } # next RCP
				
			# } # next period
			
		# } # next gcm
		
		# ### predict to FUTURE climate and FUTURE forest with deforestation OUTSIDE PAs
		# ##############################################################################
		
		# say('### predicting to future climate and future forest with deforestation outside protected areas ###')
		
		# gcmInfo <- read.csv('./Figures & Tables/Climate Model Selection/ESMs with the Greatest Difference between Current and Future Values for each BIOCLIM.csv')
		# gcms <- gcmInfo$gcm
		
		# for (gcm in gcms) {
		# # for (gcm in 'BCC-CSM1-1') {
		
			# gcmName <- gsub(gcm, pattern='-', replacement='')
			# gcmName <- tolower(gcmName)
			# gcmName <- capIt(gcmName)

			# for (period in periods) {
				
				# periodShort <- if (period == '2041 to 2060') { 2050 } else if (period == '2061 to 2080') { 2070 }

				# forestFragClassFutureName <- paste0('forestFragClass', periodShort, 'PasMasked_utm38s')
				# forestFragClassFuture <- get(forestFragClassFutureName)
				
				# for (rcp in rcps) {

					# say(gcm, ' ', period, ' ', rcp)
				
					# thisClimName <- paste0('climateFuture', gcmName, periodShort, 'Rcp', rcp)
					# futureCurrent_utm38s <- get(thisClimName)

					# preds <- stack(futureCurrent_utm38s, forestFragClassFuture)
					# prediction <- predictVareciaEnm(model, preds)
					# name <- paste0('glmEnm_vareciaGenus_climate', periodShort, gcmName, 'Rcp', rcp, '_forest', periodShort, '_defoOutsidePAs_utm38s')
					# names(prediction) <- name

					# writeRaster(prediction, paste0('./Ecological Niche Models/Prediction Rasters/', name), datatype='INT1U')
		
				# } # next RCP
				
			# } # next period
			
		# } # next gcm

# say('#####################################################################################')
# say('### create ensemble rasters for future climate ecological niche model projections ###')
# say('#####################################################################################')

	# say('Calculating mean and range of predictions across GCMs in same RCP and period.')

	# cores <- 4
	
	# ### FUTURE climate and CURRENT forest
	# #####################################
	
	# say('### ensembling future climate and current forest ###')

	# gcmInfo <- read.csv('./Figures & Tables/Climate Model Selection/ESMs with the Greatest Difference between Current and Future Values for each BIOCLIM.csv')
	# gcms <- gcmInfo$gcm
	
	# for (period in periods) {
			
		# periodShort <- if (period == '2041 to 2060') { 2050 } else if (period == '2061 to 2080') { 2070 }
			
		# for (rcp in rcps) {

			# say(period, ' ', rcp)
			
			# if (exists('predictions')) rm(predictions)
	
			# for (gcm in gcms) {
			
				# gcmName <- gsub(gcm, pattern='-', replacement='')
				# gcmName <- tolower(gcmName)
				# gcmName <- capIt(gcmName)

				# name <- paste0('glmEnm_vareciaGenus_climate', periodShort, gcmName, 'Rcp', rcp, '_forest2014_utm38s')
				# prediction <- raster(paste0('./Ecological Niche Models/Prediction Rasters/', name, '.tif'))
				
				# predictions <- if (!exists('predictions')) {
					# prediction
				# } else {
					# raster::stack(predictions, prediction)
				# }
				
			# } # next GCM

			# meanPrediction <- ensembleVareciaEnmRasters(predictions, cores=cores)
			# rm(predictions); gc()
			
			# meanName <- paste0('glmEnm_vareciaGenus_climate', periodShort, 'EnsembleMean', 'Rcp', rcp, '_forest2014_utm38s')
			# names(meanPrediction) <- meanName
			
			# writeRaster(meanPrediction, paste0('./Ecological Niche Models/Prediction Rasters/', meanName), datatype='INT1U')
			
			# rm(meanPrediction); gc()
			
		# } # next RCP
		
	# } # next period
	
	# ### FUTURE climate and FUTURE forest with deforestation ANYWHERE
	# ################################################################
	
	# say('### ensembling to future climate and future forest with deforestation anywhere ###')
	
	# gcmInfo <- read.csv('./Figures & Tables/Climate Model Selection/ESMs with the Greatest Difference between Current and Future Values for each BIOCLIM.csv')
	# gcms <- gcmInfo$gcm
	
	# for (period in periods) {
		
		# periodShort <- if (period == '2041 to 2060') { 2050 } else if (period == '2061 to 2080') { 2070 }

		# for (rcp in rcps) {

			# say(period, ' ', rcp)
			
			# if (exists('predictions')) rm(predictions)
			
			# for (gcm in gcms) {
			
				# gcmName <- gsub(gcm, pattern='-', replacement='')
				# gcmName <- tolower(gcmName)
				# gcmName <- capIt(gcmName)

				# name <- paste0('glmEnm_vareciaGenus_climate', periodShort, gcmName, 'Rcp', rcp, '_forest', periodShort, '_defoAnywhere_utm38s')
				# prediction <- raster(paste0('./Ecological Niche Models/Prediction Rasters/', name, '.tif'))
				
				# predictions <- if (!exists('predictions')) {
					# prediction
				# } else {
					# raster::stack(predictions, prediction)
				# }
	
			# } # next GCM

			# meanPrediction <- ensembleVareciaEnmRasters(predictions, cores=cores)
			
			# meanName <- paste0('glmEnm_vareciaGenus_climate', periodShort, 'EnsembleMean', 'Rcp', rcp, '_forest', periodShort, '_defoAnywhere_utm38s')
			# names(meanPrediction) <- meanName
			
			# writeRaster(meanPrediction, paste0('./Ecological Niche Models/Prediction Rasters/', meanName), datatype='INT1U')
			
			# rm(meanPrediction); gc()
			
		# } # next RCP
		
	# } # next period
	
	# ### FUTURE climate and FUTURE forest with deforestation OUTSIDE PAs
	# ###################################################################

	# say('### ensembling to future climate and future forest with deforestation outside protected areas ###')
	
	# gcmInfo <- read.csv('./Figures & Tables/Climate Model Selection/ESMs with the Greatest Difference between Current and Future Values for each BIOCLIM.csv')
	# gcms <- gcmInfo$gcm
	
	# for (period in periods) {
		
		# periodShort <- if (period == '2041 to 2060') { 2050 } else if (period == '2061 to 2080') { 2070 }

		# for (rcp in rcps) {

			# say(period, ' ', rcp)
			
			# if (exists('predictions')) rm(predictions)
			
			# for (gcm in gcms) {
			
				# gcmName <- gsub(gcm, pattern='-', replacement='')
				# gcmName <- tolower(gcmName)
				# gcmName <- capIt(gcmName)

				# name <- paste0('glmEnm_vareciaGenus_climate', periodShort, gcmName, 'Rcp', rcp, '_forest', periodShort, '_defoOutsidePAs_utm38s')
				# prediction <- raster(paste0('./Ecological Niche Models/Prediction Rasters/', name, '.tif'))
				
				# predictions <- if (!exists('predictions')) {
					# prediction
				# } else {
					# raster::stack(predictions, prediction)
				# }
	
			# } # next GCM

			# meanPrediction <- ensembleVareciaEnmRasters(predictions, cores=cores)
			
			# meanName <- paste0('glmEnm_vareciaGenus_climate', periodShort, 'EnsembleMean', 'Rcp', rcp, '_forest', periodShort, '_defoOutsidePAs_utm38s')
			# names(meanPrediction) <- meanName
			
			# writeRaster(meanPrediction, paste0('./Ecological Niche Models/Prediction Rasters/', meanName), datatype='INT1U')
			
			# rm(meanPrediction); gc()
			
		# } # next RCP
		
	# } # next period

# say('#########################################################################################')
# say('### calculate mean environmental suitability across entire region and elevation bands ###')
# say('#########################################################################################')

	# stats <- data.frame()
	# rastFiles <- listFiles('./Ecological Niche Models/Prediction Rasters', pattern='.tif')

	# # make masks of elevational bands
	# elev <- raster('./Data/Topography - GMTED2010/elevationGmted2010_utm38s.tif')
	# elevSeq <- seq(0, 100 * ceiling(maxValue(elev) / 100), by=100)
	# elevs <- list()
	# for (i in 1:(length(elevSeq) - 1)) elevs[[i]] <- list()
	# elev <- calc(elev, fun=function(x) ifelse(x < 0, 0, x))
	
	# inBand <- function(x) { ifelse(x >= bottom & x < top, 1, NA) }

	# beginCluster(4)
		
		# # by band
		# for (i in 1:(length(elevSeq) - 1)) {
			
			# bottom <- elevSeq[i]
			# top <- elevSeq[i + 1]

			# name <- paste0('elev', prefix(bottom, 4), 'to', prefix(top, 4), '_m')
			# say(name)

			# band <- clusterR(elev, calc, args=list(fun=inBand), export=c('bottom', 'top'))

			# names(band) <- name

			# # dirCreate('./Study Region & Masks/Elevational Band Masks')
			# # writeRaster(band, paste0('./Study Region & Masks/Elevational Band Masks/', name), datatype='INT1U')
			
			# elevs[[i]]$bottom <- bottom
			# elevs[[i]]$top <- top
			# elevs[[i]]$band <- band
			# elevs[[i]]$area_km2 <- cellStats(band, 'sum') * 30^2 / 1000^2

		# }

		# # save(elevs, file='./Study Region & Masks/Elevational Band Masks/Elevational Band Masks.RData')
		
	# endCluster()
	
	# # GCMs
	# gcmInfo <- read.csv('./Figures & Tables/Climate Model Selection/ESMs with the Greatest Difference between Current and Future Values for each BIOCLIM.csv')
	
	# gcms <- gcmInfo$gcm
	
	# for (rastFile in rastFiles) {
	
		# fileName <- basename(rastFile)
		# say(fileName, pre=2)
		
		# fileName <- tolower(fileName)
		
		# taxon <- if (grepl(fileName, pattern='genus')) { 'genus' } else
			# if (grepl(fileName, pattern='variegata')) { 'variegata' } else
			# if (grepl(fileName, pattern='rubra')) { 'rubra' }

		# thisGcm <- NA
		# if (grepl(fileName, pattern='current')) { thisGcm <- 'current' } else {
			# if (grepl(fileName, pattern='ensemblemean')) { thisGcm <- 'Ensemble' } else {
				
				# for (gcm in gcms) {
				
					# gcmName <- gsub(gcm, pattern='-', replacement='')
					# gcmName <- tolower(gcmName)
					# if (grepl(fileName, pattern=gcmName)) { thisGcm <- gcmName }
					
				# }
						
				# if (!is.na(thisGcm)) {
					# thisGcm <- gcms[thisGcm == tolower(gsub(gcms, pattern='-', replacement=''))]
				# }
		
			# }
			
		# }
		
		# climPeriod <- if (grepl(fileName, pattern='climatecurrent')) { 'current' } else
			# if (grepl(fileName, pattern='climate2050')) { 2050 } else
			# if (grepl(fileName, pattern='climate2070')) { 2070 }
			
		# rcp <- if (grepl(fileName, pattern='rcp2pt6')) { 2.6 } else
			# if (grepl(fileName, pattern='rcp4pt5')) { 4.5 } else
			# if (grepl(fileName, pattern='rcp6pt0')) { 6.0 } else
			# if (grepl(fileName, pattern='rcp8pt5')) { 8.5 } else { NA }
			
		# forestYear <- if (grepl(fileName, pattern='forest2014')) { 2014 } else
			# if (grepl(fileName, pattern='forest2050')) { 2050 } else
			# if (grepl(fileName, pattern='forest2070')) { 2070 } else { NA }
			
		# protection <- if (grepl(fileName, pattern='defooutsidepas')) { 'strict' } else
			# if (grepl(fileName, pattern='defoanywhere')) { 'relaxed' } else { NA }

		# # calculate mean suitability across raster
		# suitRast <- raster(rastFile)
		# meanSuit <- cellStats(suitRast, 'mean') / 100

		# # calculate mean suitability-weighted elevation
		# sumSuit <- cellStats(suitRast, 'sum') / 100
		# elevBySuit <- elev * suitRast
		# elevBySuit <- cellStats(elevBySuit, 'sum') / 100
		# meanSuitWeightedElev <- elevBySuit / sumSuit
			
		# # calculate suitability across elevational bands
		# if (exists('meanSuitByElev')) rm(meanSuitByElev)
		# for (i in seq_along(elevs)) {
		
			# suitInBand <- suitRast * elevs[[i]]$band
			# suitInBand <- cellStats(suitInBand, 'mean') / 100
			# suitInBand <- round(suitInBand, 3)
			
			# thisMeanSuitByElev <- data.frame(x = suitInBand)
			# names(thisMeanSuitByElev) <- paste0('meanSuitAtElev', prefix(elevs[[i]]$bottom, 4), 'to', prefix(elevs[[i]]$top, 4), '_m')
			
			# meanSuitByElev <- if (exists('meanSuitByElev')) {
				# cbind(meanSuitByElev, thisMeanSuitByElev)
			# } else {
				# thisMeanSuitByElev
			# }
		
			# gc()
		
		# }
			
		# thisStats <- data.frame(
			# algo = 'glm',
			# taxon = taxon,
			# climPeriod = climPeriod,
			# gcm = thisGcm,
			# rcp = rcp,
			# forestYear = forestYear,
			# protection = protection,
			# meanSuitWeightedElev_m = meanSuitWeightedElev,
			# meanSuit = round(meanSuit, 3),
			# sumSuit = round(sumSuit, 3)
		# )
		
		# thisStats <- cbind(thisStats, meanSuitByElev)
	
		# print(thisStats)

		# stats <- rbind(stats, thisStats)
			
		# rm(suitRast); gc()
			
	# } # next raster
		
	# # remember
	# dirCreate('./Figures & Tables/Ecological Niche Models - Statistics')

	# write.csv(stats, './Figures & Tables/Ecological Niche Models - Statistics/Mean Predicted Suitability by Scenario and Elevational Band.csv', row.names=FALSE)
	
	# # remember area in each elevational band
	# areaInEachElevBand_km2 <- rep(NA, length(elevs))
	# for (i in seq_along(elevs)) areaInEachElevBand_km2[i] <- elevs[[i]]$area_km2
	
	# areaInEachElevBand_km2 <- data.frame(bottom = NA, top = NA, area_km2 = areaInEachElevBand_km2)
	# for (i in seq_along(elevs)) {
		# areaInEachElevBand_km2$bottom[i] <- elevs[[i]]$bottom
		# areaInEachElevBand_km2$top[i] <- elevs[[i]]$top
	# }
	
	# write.csv(areaInEachElevBand_km2, './Figures & Tables/Ecological Niche Models - Statistics/Area in Each Elevational Band in Humid Eastern Forest Plus Buffer.csv', row.names=FALSE)
	
# say('#################################################################################################')
# say('### report values of mean environmental suitability across entire region and elevation bands ###')
# say('#################################################################################################')

	# stats <- read.csv('./Figures & Tables/Ecological Niche Models - Statistics/Mean Predicted Suitability by Scenario and Elevational Band.csv')
	
	# ### changes in suitability
	# ##########################

		# say('changes in suitability', level=2)
	
		# base <- stats$meanSuit[stats$climPeriod == 'current' & stats$gcm == 'current' & stats$forestYear == 2014 & is.na(stats$protection)]

		# ### forest loss only
		
		# # mean change in suitability assuming only forest loss with STRICT protection and NO CLIMATE CHANGE by 2070
		# vals <- stats$meanSuit[stats$climPeriod == 'current' & stats$gcm == 'current' & stats$forestYear == 2070 & stats$protection == 'strict']
		# delta <- -100 * (base - vals) / base
		# say('Change in mean suitability due only to forest loss by 2070 assuming STRICT protection: ', sprintf('%.0f', delta), '%')

		# # mean change in suitability assuming only forest loss with RELAXED protection and NO CLIMATE CHANGE by 2070
		# vals <- stats$meanSuit[stats$climPeriod == 'current' & stats$gcm == 'current' & stats$forestYear == 2070 & stats$protection == 'relaxed']
		# delta <- -100 * (base - vals) / base
		# say('Change in mean suitability due only to forest loss by 2070 assuming RELAXED protection: ', sprintf('%.0f', delta), '%')

		# ### climate change only
		
		# # mean change in suitability assuming NO FOREST LOSS with but with CLIMATE CHANGE by 2070
		# vals <- stats$meanSuit[stats$climPeriod == '2070' & stats$gcm != 'current' & stats$gcm != 'Ensemble' & stats$rcp == 8.5 & stats$forestYear == 2014]
		# delta <- -100 * mean((base - vals) / base)
		# lower <- -100 * (base - min(vals)) / base
		# upper <- -100 * (base - max(vals)) / base
		# say('Change in mean suitability due only to climate change under RCP8.5 by 2070 assuming NO FOREST LOSS: ', sprintf('%.0f', delta), '% (range: ', sprintf('%.0f', upper), ' to ', sprintf('%.0f', lower), '%)')

		# ### climate change and forest loss
		
		# # mean change in suitability assuming FOREST LOSS WITH STRICT protection and CLIMATE CHANGE by 2070
		# vals <- stats$meanSuit[stats$climPeriod == '2070' & stats$gcm != 'current' & stats$gcm != 'Ensemble' & stats$rcp == 8.5 & stats$forestYear == 2070 & stats$protection == 'strict']
		# delta <- -100 * mean((base - vals) / base)
		# lower <- -100 * (base - min(vals)) / base
		# upper <- -100 * (base - max(vals)) / base
		# say('Change in mean suitability due to FOREST LOSS assuming STRICT protection and CLIMATE CHANGE under RCP8.5: ', sprintf('%.0f', delta), '% (range: ', sprintf('%.0f', upper), ' to ', sprintf('%.0f', lower), '%)')

		# # mean change in suitability assuming FOREST LOSS WITH RELAXED protection and CLIMATE CHANGE by 2070
		# vals <- stats$meanSuit[stats$climPeriod == '2070' & stats$gcm != 'current' & stats$gcm != 'Ensemble' & stats$rcp == 8.5 & stats$forestYear == 2070 & stats$protection == 'relaxed']
		# delta <- -100 * mean((base - vals) / base)
		# lower <- -100 * (base - min(vals)) / base
		# upper <- -100 * (base - max(vals)) / base
		# say('Change in mean suitability due to FOREST LOSS assuming RELAXED protection and CLIMATE CHANGE under RCP8.5: ', sprintf('%.0f', delta), '% (range: ', sprintf('%.0f', upper), ' to ', sprintf('%.0f', lower), '%)')

	# ### changes in elevation
	# ########################

		# say('changes in elevation', level=2)
	
		# base <- stats$meanSuitWeightedElev_m[stats$climPeriod == 'current' & stats$gcm == 'current' & stats$forestYear == 2014 & is.na(stats$protection)]

		# ### forest loss only
		
		# # mean change in elevation assuming only forest loss with STRICT protection and NO CLIMATE CHANGE by 2070
		# vals <- stats$meanSuitWeightedElev_m[stats$climPeriod == 'current' & stats$gcm == 'current' & stats$forestYear == 2070 & stats$protection == 'strict']
		# delta <- vals - base
		# say('Change in mean elevation due only to forest loss by 2070 assuming STRICT protection: ', sprintf('%.0f', delta), ' m')

		# # mean change in elevation assuming only forest loss with RELAXED protection and NO CLIMATE CHANGE by 2070
		# vals <- stats$meanSuitWeightedElev_m[stats$climPeriod == 'current' & stats$gcm == 'current' & stats$forestYear == 2070 & stats$protection == 'relaxed']
		# delta <- vals - base
		# say('Change in mean elevation due only to forest loss by 2070 assuming RELAXED protection: ', sprintf('%.0f', delta), ' m')

		# ### climate change only
		
		# # mean change in elevation assuming NO FOREST LOSS with but with CLIMATE CHANGE by 2070
		# vals <- stats$meanSuitWeightedElev_m[stats$climPeriod == '2070' & stats$gcm != 'current' & stats$rcp == 8.5 & stats$forestYear == 2014]
		# delta <- mean(vals - base)
		# say('Change in mean elevation due only to climate change under RCP8.5 by 2070 assuming NO FOREST LOSS: ', sprintf('%.0f', delta), ' m')

		# ### climate change and forest loss
		
		# # mean change in elevation assuming FOREST LOSS WITH STRICT protection and CLIMATE CHANGE by 2070
		# vals <- stats$meanSuitWeightedElev_m[stats$climPeriod == '2070' & stats$gcm != 'current' & stats$rcp == 8.5 & stats$forestYear == 2070 & stats$protection == 'strict']
		# delta <- mean(vals - base)
		# say('Change in mean elevation due to FOREST LOSS assuming STRICT protection and CLIMATE CHANGE under RCP8.5: ', sprintf('%.0f', delta), ' m')

		# # mean change in elevation assuming FOREST LOSS WITH RELAXED protection and CLIMATE CHANGE by 2070
		# vals <- stats$meanSuitWeightedElev_m[stats$climPeriod == '2070' & stats$gcm != 'current' & stats$rcp == 8.5 & stats$forestYear == 2070 & stats$protection == 'relaxed']
		# delta <- mean(vals - base)
		# say('Change in mean elevation due to FOREST LOSS assuming RELAXED protection and CLIMATE CHANGE under RCP8.5: ', sprintf('%.0f', delta), ' m')

# say('##################################################################')
# say('### compare elevational distribution of forest and occurrences ###')
# say('##################################################################')

	# say('Wanting to determine if current distribution of forest seems to curtail apparent climatic niche.')

	# humidForestBufferMask_utm38s <- raster('./Study Region & Masks/UTM 38S 30-m Resolution/humidForestBufferMask_utm38s.tif')
	
	# # current forest
	# forest2014 <- raster('./Data/Forest - Vieilledent et al 2018/forest2014.tif')
	# forest2014 <- crop(forest2014, humidForestBufferMask_utm38s)
	# gc()
	
	# # current forest fragmentation
	# frag2014 <- raster('./Data/Forest - Vieilledent et al 2018/forestFragClass2014_utm38s.tif')
	# frag2014 <- crop(frag2014, humidForestBufferMask_utm38s)
	# interior2014 <- (frag2014 == 6)
	# gc()
	
	# # elevation
	# elev <- raster('./Data/Topography - GMTED2010/elevationGmted2010_utm38s.tif')
	# elevSeq <- seq(0, 100 * ceiling(maxValue(elev) / 100), by=100)
	
	# inBand <- function(x) { ifelse(x >= bottom & x < top, 1, NA) }

	# # survey sites
	# load('./Ecological Niche Models/Collated Presence and Background Point Data for Varecia GENUS.RData')
	# v <- taxonData[taxonData$presBg == 1, ]
	# load('./Ecological Niche Models/Collated Presence and Background Point Data for Varecia VARIEGATA.RData')
	# vv <- taxonData[taxonData$presBg == 1, ]
	# load('./Ecological Niche Models/Collated Presence and Background Point Data for Varecia RUBRA.RData')
	# vr <- taxonData[taxonData$presBg == 1, ]

	# vSp <- SpatialPoints(v[ , longLat], getCRS('wgs84', TRUE))
	# vvSp <- SpatialPoints(vv[ , longLat], getCRS('wgs84', TRUE))
	# vrSp <- SpatialPoints(vr[ , longLat], getCRS('wgs84', TRUE))
	
	# vSp <- sp::spTransform(vSp, CRS(madEaProj))
	# vvSp <- sp::spTransform(vvSp, CRS(madEaProj))
	# vrSp <- sp::spTransform(vrSp, CRS(madEaProj))
	
	# v$elevation_utm38s <- raster::extract(elev, vSp)
	# vv$elevation_utm38s <- raster::extract(elev, vvSp)
	# vr$elevation_utm38s <- raster::extract(elev, vrSp)
	
	# # for storing proportion of survey sites, elevation, and forest cover in each band
	# stats <- data.frame()
	
	# # cores and excluded packages (not needed on nodes)
	# cores <- 3
	# exclude <- c('omnibus', 'enmSdm', 'statisfactory', 'legendary', 'fasterRaster', 'dismo', 'rgeos', 'geosphere', 'brglm2', 'phcfM', 'fpCompare', 'scales', 'rgl', 'rayshader', 'tictoc')
	
	# # convert elevations <0 to 0
	# lt0to0 <- function(x) ifelse(x < 0, 0, x)
	# beginCluster(cores, type='SOCK', exclude=exclude)
		# elev <- clusterR(elev, calc, args=list(fun=lt0to0))
	# endCluster()
	
		# # tally area, forest area, and occurrences by elevational band
		# for (i in 1:(length(elevSeq) - 1)) {
			
			# bottom <- elevSeq[i]
			# top <- elevSeq[i + 1]

			# name <- paste0('elev', prefix(bottom, 4), 'to', prefix(top, 4), '_m')
			# say(name)

			# beginCluster(cores, type='SOCK', exclude=exclude)
				# elevBand <- clusterR(elev, calc, args=list(fun=inBand), export=c('bottom', 'top'))
				# gc()
			# endCluster()
			# gc()
			
			# forestBand <- forest2014 * elevBand
			# gc()
			
			# interiorBand <- interior2014 * elevBand
			# gc()

			# elevArea_km2 <- cellStats(elevBand, 'sum') * 30^2 / 1000^2
			# forestArea_km2 <- cellStats(forestBand, 'sum') * 30^2 / 1000^2
			# forestInteriorArea_km2 <- cellStats(interiorBand, 'sum') * 30^2 / 1000^2

			# rm(elevBand, forestBand, interiorBand); gc()
			
			# genusInBand <- sum(v$elevation_utm38s >= bottom & v$elevation_utm38s < top)
			# vvInBand <- sum(vv$elevation_utm38s >= bottom & vv$elevation_utm38s < top)
			# vrInBand <- sum(vr$elevation_utm38s >= bottom & vr$elevation_utm38s < top)
			
			# thisStats <- data.frame(
				# bottom = bottom,
				# top = top,
				# elevArea_km2 = elevArea_km2,
				# forestArea_km2 = forestArea_km2,
				# forestInteriorArea_km2 = forestInteriorArea_km2,
				# genusInBand = genusInBand,
				# vvInBand = vvInBand,
				# vrInBand = vrInBand
			# )
			
			# stats <- rbind(stats, thisStats)
			
		# }

	# endCluster()

	# rownames(stats) <- 1:nrow(stats)
	
	# dirCreate('./Figures & Tables/Occurrences vs Forest Cover by Elevation')
	# write.csv(stats, './Figures & Tables/Occurrences vs Forest Cover by Elevation/Occurrences vs Forest Cover by Elevation.csv')

	# ### plot
	# stats <- read.csv('./Figures & Tables/Occurrences vs Forest Cover by Elevation/Occurrences vs Forest Cover by Elevation.csv')
	
	# elev <- raster('./Data/Topography - GMTED2010/elevationGmted2010_utm38s.tif')
	# elevSeq <- seq(0, 100 * ceiling(maxValue(elev) / 100), by=100)

	# elevCumSum <- cumsum(stats$elevArea_km2)
	# forestCumSum <- cumsum(stats$forestArea_km2)
	# interiorCumSum <- cumsum(stats$forestInteriorArea_km2)
	# genusInBandCumSum <- cumsum(stats$genusInBand)
	# vvInBandCumSum <- cumsum(stats$vvInBand)
	# vrInBandCumSum <- cumsum(stats$vrInBand)
	
	# elevCumSum <- elevCumSum / max(elevCumSum)
	# forestCumSum <- forestCumSum / max(forestCumSum)
	# interiorCumSum <- interiorCumSum / max(interiorCumSum)
	# genusInBandCumSum <- genusInBandCumSum / max(genusInBandCumSum)
	# vvInBandCumSum <- vvInBandCumSum / max(vvInBandCumSum)
	# vrInBandCumSum <- vrInBandCumSum / max(vrInBandCumSum)
	
	# x <- rowMeans(stats[ , c('bottom', 'top')])

	# png('./Figures & Tables/Cumulative Occupancy of Available Habitat.png', width=800, height=800, res=200)
		
		# par(oma=0.2 * c(1, 1, 1, 1), mar=0.5 * c(8, 8, 4, 2), cex.axis=0.8)
		
		# plot(x, elevCumSum, ylim=c(0, 1), type='l', lwd=2, xlab='Elevation (m)', ylab='Cumulative proportion')
		# lines(x, forestCumSum, lwd=2, col='forestgreen')
		# lines(x, interiorCumSum, lwd=2, lty='dashed', col='forestgreen')
		# lines(x, genusInBandCumSum, lwd=3, col='red')
		# lines(x, vvInBandCumSum, lwd=2, col='blue')
		# lines(x, vrInBandCumSum, lwd=2, col='orange')

		# legend('bottomright', inset=0.01, legend=c('Area', 'Forest cover', 'Interior forest', 'Genus occurrences', 'V. variegata occurrences', 'V. rubra occurrences'), col=c('black', 'forestgreen', 'forestgreen', 'red', 'blue', 'orange'), lwd=2, lty=c('solid', 'solid', 'dashed', 'solid', 'solid', 'solid'), cex=0.6)
		
	# dev.off()
	
# say('##########################################################')
# say('### evaluate changes in suitability in protected areas ###')
# say('##########################################################')

	# say('Wanting to characterize current and potential future suitability to Varecia in each protected area.')

	# # load ENM raster
	# sq <- raster('./Ecological Niche Models/Prediction Rasters/glmEnm_vareciaGenus_climateCurrent_forest2014_utm38s.tif')

	# fut_2050rcp26_strict <- raster('./Ecological Niche Models/Prediction Rasters/glmEnm_vareciaGenus_climate2050EnsembleMeanRcp2pt6_forest2050_defoOutsidePAs_utm38s.tif')
	# fut_2050rcp26_relaxed <- raster('./Ecological Niche Models/Prediction Rasters/glmEnm_vareciaGenus_climate2050EnsembleMeanRcp2pt6_forest2050_defoAnywhere_utm38s.tif')
	
	# fut_2070rcp26_strict <- raster('./Ecological Niche Models/Prediction Rasters/glmEnm_vareciaGenus_climate2070EnsembleMeanRcp2pt6_forest2070_defoOutsidePAs_utm38s.tif')
	# fut_2070rcp26_relaxed <- raster('./Ecological Niche Models/Prediction Rasters/glmEnm_vareciaGenus_climate2070EnsembleMeanRcp2pt6_forest2070_defoAnywhere_utm38s.tif')
	
	# fut_2050rcp45_strict <- raster('./Ecological Niche Models/Prediction Rasters/glmEnm_vareciaGenus_climate2050EnsembleMeanRcp4pt5_forest2050_defoOutsidePAs_utm38s.tif')
	# fut_2050rcp45_relaxed <- raster('./Ecological Niche Models/Prediction Rasters/glmEnm_vareciaGenus_climate2050EnsembleMeanRcp4pt5_forest2050_defoAnywhere_utm38s.tif')
	
	# fut_2070rcp45_strict <- raster('./Ecological Niche Models/Prediction Rasters/glmEnm_vareciaGenus_climate2070EnsembleMeanRcp4pt5_forest2070_defoOutsidePAs_utm38s.tif')
	# fut_2070rcp45_relaxed <- raster('./Ecological Niche Models/Prediction Rasters/glmEnm_vareciaGenus_climate2070EnsembleMeanRcp4pt5_forest2070_defoAnywhere_utm38s.tif')
	
	# fut_2050rcp60_strict <- raster('./Ecological Niche Models/Prediction Rasters/glmEnm_vareciaGenus_climate2050EnsembleMeanRcp6pt0_forest2050_defoOutsidePAs_utm38s.tif')
	# fut_2050rcp60_relaxed <- raster('./Ecological Niche Models/Prediction Rasters/glmEnm_vareciaGenus_climate2050EnsembleMeanRcp6pt0_forest2050_defoAnywhere_utm38s.tif')
	
	# fut_2070rcp60_strict <- raster('./Ecological Niche Models/Prediction Rasters/glmEnm_vareciaGenus_climate2070EnsembleMeanRcp6pt0_forest2070_defoOutsidePAs_utm38s.tif')
	# fut_2070rcp60_relaxed <- raster('./Ecological Niche Models/Prediction Rasters/glmEnm_vareciaGenus_climate2070EnsembleMeanRcp6pt0_forest2070_defoAnywhere_utm38s.tif')
	
	# fut_2050rcp85_strict <- raster('./Ecological Niche Models/Prediction Rasters/glmEnm_vareciaGenus_climate2050EnsembleMeanRcp8pt5_forest2050_defoOutsidePAs_utm38s.tif')
	# fut_2050rcp85_relaxed <- raster('./Ecological Niche Models/Prediction Rasters/glmEnm_vareciaGenus_climate2050EnsembleMeanRcp8pt5_forest2050_defoAnywhere_utm38s.tif')
	
	# fut_2070rcp85_strict <- raster('./Ecological Niche Models/Prediction Rasters/glmEnm_vareciaGenus_climate2070EnsembleMeanRcp8pt5_forest2070_defoOutsidePAs_utm38s.tif')
	# fut_2070rcp85_relaxed <- raster('./Ecological Niche Models/Prediction Rasters/glmEnm_vareciaGenus_climate2070EnsembleMeanRcp8pt5_forest2070_defoAnywhere_utm38s.tif')
	
	# futs <- c('fut_2050rcp26_strict', 'fut_2070rcp26_strict', 'fut_2050rcp26_relaxed', 'fut_2070rcp26_relaxed', 'fut_2050rcp45_strict', 'fut_2070rcp45_strict', 'fut_2050rcp45_relaxed', 'fut_2070rcp45_relaxed', 'fut_2050rcp60_strict', 'fut_2070rcp60_strict', 'fut_2050rcp60_relaxed', 'fut_2070rcp60_relaxed', 'fut_2050rcp85_strict', 'fut_2070rcp85_strict', 'fut_2050rcp85_relaxed', 'fut_2070rcp85_relaxed')
	
	# # eastern humid forest
	# humidForestBufferMask_utm38s <- raster('./Study Region & Masks/UTM 38S 30-m Resolution/humidForestBufferMask_utm38s.tif')
	
	# # protected area
	# pas <- shapefile('./Data/Protected Areas/WDPA_Sept2018_MDG-shapefile-polygons')
	# pas <- sp::spTransform(pas, CRS(madEaProj))
	
	# pas <- crop(pas, humidForestBufferMask_utm38s)
	# pas$area_km2 <- gArea(pas, byid=TRUE) / 1000^2
	
	# paNames <- sort(unique(pas$NAME))
	
	# # for output
	# stats <- data.frame()
	
	# # by PA
	# for (thisPa in paNames) {
	
		# say(thisPa)
	
		# thisPaPoly <- pas[pas$NAME == thisPa, ]
	
		# # mask
		# paMask <- crop(humidForestBufferMask_utm38s, thisPaPoly)
		# paMask <- rasterize(thisPaPoly, paMask)
	
		# # area of PA in eastern humid forest + buffer
		# areaInHumidForestBuffer_km2 <- cellStats(paMask, 'sum') * 30^2 / 1000^2
	
		# # present suitability
		# paSq <- crop(sq, thisPaPoly)
		# paSq <- paMask * paSq
		# suitSq <- cellStats(paSq, 'sum') / 100
		
		# # remember
		# thisStats <- data.frame(
			# pa = thisPa,
			# iucnCategory = thisPaPoly$IUCN_CAT,
			# area_km2 = thisPaPoly$area_km2,
			# areaInHumidForestBuffer_km2 = areaInHumidForestBuffer_km2,
			# presentSuit_sum = suitSq
		# )
		
		# # future suitability
		# for (thisFut in futs) {
		
			# period <- substr(thisFut, 5, 8)
			# rcp <- substr(thisFut, 12, 13)
			# prot <- substr(thisFut, 15, nchar(thisFut))
		
			# x <- get(thisFut)
			
			# paFut <- crop(x, thisPaPoly)
			# paFut <- paMask * paFut
			# suitFut <- cellStats(paFut, 'sum') / 100
		
			# thisOut <- data.frame(x=suitFut)
			# names(thisOut) <- thisFut

			# thisStats <- cbind(thisStats, thisOut)
			
		# }
		
		# stats <- rbind(stats, thisStats)
		
	# }
	
	# write.csv(stats, './Figures & Tables/Ecological Niche Models - Statistics/Change in Suitability by Protected Area.csv', row.names=FALSE)
	
# say('##########################################')
# say('### create fine-scale hillshade raster ###')
# say('##########################################')

	# ### all of Madagascar
		
	# load('./Study Region & Masks/UTM 38S 30-m Resolution/Eastern Humid Forest Polygon.RData')
	# load('./Study Region & Masks/UTM 38S 30-m Resolution/Madagascar from GADM 3.6.RData')

	# elev <- raster('./Data/Topography - GMTED2010/elevationGmted2010_wgs84.tif')

	# madagascarBuffer10km_utm38s <- gBuffer(madagascar_utm38s, width=10000)
	# madagascar_wgs84 <- sp::spTransform(madagascar_utm38s, getCRS('wgs84', TRUE))
	
	# elev <- crop(elev, madagascar_wgs84)

	# beginCluster(3)
	# elev <- projectRaster(elev, crs=madEaProj)
	# endCluster()
	
	# aspect <- terrain(elev, 'aspect')
	# slope <- terrain(elev, 'slope')
	
	# names(aspect) <- 'aspectGmted2010_utm38s'
	# names(slope) <- 'slopeGmted2010_utm38s'

	# hs <- hillShade(slope=slope, aspect=aspect, direction=60, normalize=TRUE)
	
	# mask <- rasterize(madagascar_utm38s, hs)
	# mask <- 1 + 0 * mask
	# hs <- hs * mask
	# hs <- round(hs)

	# names(hs) <- 'hillshadeGmted2010_utm38s'
	# writeRaster(hs, './Data/Topography - GMTED2010/hillshadeGmted2010_utm38s', datatype='INT1U')
	
# say('#############################################################')
# say('### create displays of ecological niche model predictions ###')
# say('#############################################################')

	# # generalization
	# gcm <- 'EnsembleMean'
		
	# outDir <- './Figures & Tables/Ecological Niche Models - Predictions/'
	# dirCreate(outDir)

	# # cexMult <- 0.45 # multiplier for cex... good for PNG
	# cexMult <- 4 * 0.45 # multiplier for cex

	# # RCPs for each column
	# rcpCol2 <- '2pt6'
	# rcpCol3 <- '4pt5'
	# rcpCol4 <- '6pt0'
	# rcpCol5 <- '8pt5'
	
	# rcpNice <- function(rcp) gsub(rcp, pattern='pt', replacement='.')
	# rcpCol2Nice <- rcpNice(rcpCol2)
	# rcpCol3Nice <- rcpNice(rcpCol3)
	# rcpCol4Nice <- rcpNice(rcpCol4)
	# rcpCol5Nice <- rcpNice(rcpCol5)
		
	# # ancillary geo data
	# humidForestBufferMask_utm38s <- raster('./Study Region & Masks/UTM 38S 30-m Resolution/humidForestBufferMask_utm38s.tif')

	# load('./Study Region & Masks/UTM 38S 30-m Resolution/Eastern Humid Forest Polygon.RData')
	# load('./Study Region & Masks/UTM 38S 30-m Resolution/Madagascar from GADM 3.6.RData')
	
	# # # pas <- shapefile('./Data/Protected Areas/WDPA_Sept2018_MDG-shapefile-polygons')
	# # # pas <- sp::spTransform(pas, CRS(madEaProj))
	# # # pas <- crop(pas, madagascar_utm38s)
	
	# # load('./Data/Protected Areas/WDPA_Sept2018_MDG-shapefile-polygons-onlyTerrestrial.RData')
	# load('./Data/Protected Areas/WDPA_Sept2018_MDG-shapefile-polygons_easternRainforest.RData')
	# pas_easternRainforest <- sp::spTransform(pas_easternRainforest, CRS(madEaProj))
	
	# # list of maps
	# maps <- c('r1c1', 'r2c1', 'r3c1', 'r1c2', 'r2c2', 'r3c2', 'r1c3', 'r2c3', 'r3c3', 'r1c4', 'r2c4', 'r3c4', 'r1c5', 'r2c5', 'r3c5')
	# # maps <- c('r1c1', 'r2c1')
	# # maps <- c('r3c5')
	
	# # color ramp for focal rasters
	# cols <- colorRampPalette(c('lightsalmon', 'red', 'darkred'))
	# colBreaks <- 5 # number of color categories plus 1
	# cols <- cols(colBreaks - 2)
	# cols <- c(NA, cols)

	# colsLegend <- cols
	# colsLegend[is.na(colsLegend)] <- 'white'

	# # hillshade
	# hs <- raster('./Data/Topography - GMTED2010/hillshadeGmted2010_utm38s.tif')
	
	# # hillshade colors
	# hsCols <- paste0('gray', 0:100)
	# hsCols <- alpha(hsCols, 0.7)

	# ### setup graticule
	# longs <- seq(45, 51, by=2)
	# lats <- seq(-12, -26, by=-2)
	
	# longLim <- range(longs)# + c(-0.02, 0.02)
	# latLim <- range(lats)# + c(-0.02, 0.02)
	
	# grats <- graticule(longs, lats, proj=madEaProj, xlim=longLim, ylim=latLim)
	
	# for (futYear in c(2050, 2070)) {
	# # for (futYear in c(2070)) {
	
		# say(futYear)

		# # prediction rasters
		# # rasters (u = upper, m = middle, l = lower or left, r = right)
		# rastDir <- paste0('./Ecological Niche Models/Prediction Rasters')
		
		# r1c1 <- raster(paste0(rastDir, '/glmEnm_vareciaGenus_climateCurrent_forest2014_utm38s.tif'))
		# r2c1 <- raster(paste0(rastDir, '/glmEnm_vareciaGenus_climateCurrent_forest', futYear, '_defoOutsidePAs_utm38s.tif'))
		# r3c1 <- raster(paste0(rastDir, '/glmEnm_vareciaGenus_climateCurrent_forest', futYear, '_defoAnywhere_utm38s.tif'))

		# r1c2 <- raster(paste0(rastDir, '/glmEnm_vareciaGenus_climate', futYear, gcm, 'Rcp', rcpCol2, '_forest2014_utm38s.tif'))
		# r2c2 <- raster(paste0(rastDir, '/glmEnm_vareciaGenus_climate', futYear, gcm, 'Rcp', rcpCol2, '_forest', futYear, '_defoOutsidePAs_utm38s.tif'))
		# r3c2 <- raster(paste0(rastDir, '/glmEnm_vareciaGenus_climate', futYear, gcm, 'Rcp', rcpCol2, '_forest', futYear, '_defoAnywhere_utm38s.tif'))

		# r1c3 <- raster(paste0(rastDir, '/glmEnm_vareciaGenus_climate', futYear, gcm, 'Rcp', rcpCol3, '_forest2014_utm38s.tif'))
		# r2c3 <- raster(paste0(rastDir, '/glmEnm_vareciaGenus_climate', futYear, gcm, 'Rcp', rcpCol3, '_forest', futYear, '_defoOutsidePAs_utm38s.tif'))
		# r3c3 <- raster(paste0(rastDir, '/glmEnm_vareciaGenus_climate', futYear, gcm, 'Rcp', rcpCol3, '_forest', futYear, '_defoAnywhere_utm38s.tif'))

		# r1c4 <- raster(paste0(rastDir, '/glmEnm_vareciaGenus_climate', futYear, gcm, 'Rcp', rcpCol4, '_forest2014_utm38s.tif'))
		# r2c4 <- raster(paste0(rastDir, '/glmEnm_vareciaGenus_climate', futYear, gcm, 'Rcp', rcpCol4, '_forest', futYear, '_defoOutsidePAs_utm38s.tif'))
		# r3c4 <- raster(paste0(rastDir, '/glmEnm_vareciaGenus_climate', futYear, gcm, 'Rcp', rcpCol4, '_forest', futYear, '_defoAnywhere_utm38s.tif'))

		# r1c5 <- raster(paste0(rastDir, '/glmEnm_vareciaGenus_climate', futYear, gcm, 'Rcp', rcpCol5, '_forest2014_utm38s.tif'))
		# r2c5 <- raster(paste0(rastDir, '/glmEnm_vareciaGenus_climate', futYear, gcm, 'Rcp', rcpCol5, '_forest', futYear, '_defoOutsidePAs_utm38s.tif'))
		# r3c5 <- raster(paste0(rastDir, '/glmEnm_vareciaGenus_climate', futYear, gcm, 'Rcp', rcpCol5, '_forest', futYear, '_defoAnywhere_utm38s.tif'))

		# pdf(paste0(outDir, '/ENM Predictions Using Ensemble Mean across GCMs for ', futYear, ' x Deforestation x RCPs.pdf'), width=2 * 4.40, height=1.5 * 9.00)
		
			# pars <- par(mfrow=c(3, 5), oma=c(2, 4.2, 6, 3), mar=1.2 * c(0, 0.2, 0.2, 0))

			# for (map in maps) {
			# # for (map in maps[1]) {
			# # for (map in 'r3c5') {
				
				# say('subplot ', map)
				
				# rowCol <- if (map == 'r1c1') { c(1, 1)} else
					# if (map == 'r1c2') { c(1, 2)} else
					# if (map == 'r1c3') { c(1, 3)} else
					# if (map == 'r1c4') { c(1, 4)} else
					# if (map == 'r1c5') { c(1, 5)} else
				
					# if (map == 'r2c1') { c(2, 1)} else
					# if (map == 'r2c2') { c(2, 2)} else
					# if (map == 'r2c3') { c(2, 3)} else
					# if (map == 'r2c4') { c(2, 4)} else
					# if (map == 'r2c5') { c(2, 5)} else
				
					# if (map == 'r3c1') { c(3, 1)} else
					# if (map == 'r3c2') { c(3, 2)} else
					# if (map == 'r3c3') { c(3, 3)} else
					# if (map == 'r3c4') { c(3, 4)} else
					# if (map == 'r3c5') { c(3, 5)}
				
				# letter <- if (is.na(map)) { NA } else
					# if (map == 'r1c1') { 'a' } else
					# if (map == 'r1c2') { 'b' } else
					# if (map == 'r1c3') { 'c' } else
					# if (map == 'r1c4') { 'd' } else
					# if (map == 'r1c5') { 'e' } else

					# if (map == 'r2c1') { 'f' } else
					# if (map == 'r2c2') { 'g' } else
					# if (map == 'r2c3') { 'h' } else
					# if (map == 'r2c4') { 'i' } else
					# if (map == 'r2c5') { 'j' } else
				
					# if (map == 'r3c1') { 'k' } else
					# if (map == 'r3c2') { 'l' } else
					# if (map == 'r3c3') { 'm' } else
					# if (map == 'r3c4') { 'n' } else
					# if (map == 'r3c5') { 'o' }
				
				# par(mfg=rowCol)
				# plot.new()
				
				# x <- get(map)
				
				# plot(humidForest_utm38s, col=NA, border='black', lty='solid', lwd=0.2)

				# ### graticule
				# plot(grats, add=TRUE, lwd=0.2, col='gray')

				# ### graticule longitude labels
				# if (grepl(map, pattern='r3')) {
					
					# longs <- seq(47, 51, by=2)
					# lats <- rep(-25.5, length(longs))
					
					# longLim <- range(longs)
					# latLim <- range(lats)
					
					# gratLabs <- graticule_labels(longs, lats, xline=min(longLim), yline=max(latLim), proj=madEaProj)
					# longLabs <- gratLabs[gratLabs$islon, ]
					# text(longLabs, label=parse(text=longLabs$lab), col='gray', cex=cexMult * 0.6, srt=90, xpd=NA, adj=c(0.9, 0.5))
					
				# }

				# ### main plot
				# plot(madagascar_utm38s, lwd=0.4, col='white', add=TRUE)
				# # plot(pas_easternRainforest, border=NA, col='lightskyblue', add=TRUE)
				# plot(pas_easternRainforest, border=NA, col=alpha('lightskyblue', 0.8), lwd=0.3, add=TRUE)
				# plot(x, col=cols, breaks=seq(0, 100, length.out=colBreaks), legend=FALSE, add=TRUE)
				# plot(humidForest_utm38s, border='black', lty='solid', lwd=0.2, add=TRUE)
							
				# ### graticule latitude labels
				# if (grepl(map, pattern='c5')) {
					
					# lats <- seq(-14, -24, by=-2)
					# longs <- rep(51, length(lats))
					
					# longLim <- range(longs)
					# latLim <- range(lats)
					
					# gratLabs <- graticule_labels(longs, lats, xline=min(longLim), yline=max(latLim), proj=madEaProj)
					# latLabs <- gratLabs[!gratLabs$islon, ]
					
					# labs <- parse(text=latLabs$lab)
					# text(latLabs, label=labs, col='gray', cex=cexMult * 0.6, xpd=NA, adj=c(-0.13, 0.5))
					
				# }

				# ### title
				# usr <- par('usr')
				# titleX <- usr[1] - 0 * (usr[2] - usr[1])
				# titleY <- usr[4] - 0.05 * (usr[4] - usr[3])
				# text(titleX, titleY, labels=letter, pos=4, cex=cexMult * 1.1, xpd=NA, font=2)

				# # legend
				# if (map == 'r3c5') {

					# # legend
					# labs <- c(0.25, 0.5, 0.75)
					# paSwatch <- list(swatchAdjY=c(0, 0.13), col='lightskyblue', border='black', labels='PA', lwd=0.2)
					
					# legendBreaks('bottomright', inset=c(0.17, 0.005), width=0.2, height=0.43, labels=labs, labAdjX=1.2, labAdjY=labs, cex=cexMult * 0.72, adjX=c(0, 0.7), adjY=c(0.22, 0.75), col=colsLegend, colBorder=NA, title='', boxBorder=NA, boxBg=NA, xpd=NA, swatches=list(paSwatch), lwd=0.2)
					
					# suitX <- 855000
					# suitY <- 7412000
					
					# text(suitX, suitY, srt=90, labels='   Suitability', cex=cexMult * 0.8)
					
				# }

				# ### insets
				
				# insetLwd <- 0.35
				# plot(madFocus1, lwd=insetLwd, lend=1, add=TRUE)
				# plot(madFocus2, lwd=insetLwd, lend=1, add=TRUE)
				# plot(madFocus3, lwd=insetLwd, lend=1, add=TRUE)

				# # column label (RCPs)
				# if (map %in% c('r1c1', 'r1c2', 'r1c3', 'r1c4', 'r1c5')) {

					# labCol <- if (is.na(map)) {
						# NA
					# } else if (map == 'r1c1') {
						# 'Current climate'
					# } else if (map == 'r1c2') {
						# paste0(futYear, ' climate\nRCP ', rcpCol2Nice)
					# } else if (map == 'r1c3') {
						# paste0(futYear, ' climate\nRCP ', rcpCol3Nice)
					# } else if (map == 'r1c4') {
						# paste0(futYear, ' climate\nRCP ', rcpCol4Nice)
					# } else if (map == 'r1c5') {
						# paste0(futYear, ' climate\nRCP ', rcpCol5Nice)
					# }
					
					# x <- 0.5 * (usr[1] + usr[2])
					# y <- usr[4] + 0.1 * (usr[4] - usr[3])
					# text(x, y, labels=labCol, cex=cexMult * 1, xpd=NA)
				
				# }
				
				# # row labels (defo)
				# if (map %in% c('r1c1', 'r2c1', 'r3c1')) {

				# usr <- par('usr')
					# labRow <- if (map == 'r1c1') {
						# '2014 forest'
					# } else if (map == 'r2c1') {
						# paste0(futYear, ' forest\nstrict protection')
					# } else if (map == 'r3c1') {
						# paste0(futYear, ' forest\nrelaxed protection')
					# }
				
					# x <- usr[1] - 0.22 * (usr[2] - usr[1])
					# y <- 0.5 * (usr[3] + usr[4])
					# text(x, y, labels=labRow, cex=cexMult * 1.1, xpd=NA, srt=90)
					
				# }
					
			# } # next map
			
		# dev.off()

	# } # next future year

	# # # ### insets
	# # # ##########

		# # # for (inset in 1:3) {
		# # # # for (inset in 1) {

			# # # say('inset ', inset, level=3)
		
			# # # # get bounding box
			# # # thisFocus <- paste0('madFocus', inset)

			# # # focus <- get(thisFocus)
			# # # ext <- extent(focus)
			# # # ratio <- (ext@ymax - ext@ymin)/ (ext@xmax - ext@xmin)
	
			# # # # crop geo data
			# # # humidForestBufferMaskCrop_utm38s <- crop(humidForestBufferMask_utm38s, focus)
			# # # madagascarCrop_utm38s <- crop(madagascar_utm38s, focus)
			# # # pasCrop <- crop(pas, focus)
			# # # hsCrop <- crop(hs, focus)
						
			# # # height <- 3600
			# # # width <- round(height / ratio)

			# # # png(paste0(outDir, '/ENM Predictions Using Ensemble Mean across GCMs for ', futYear, ' x Deforestation x RCPs Inset ', rcpMiddleRow, ' & ', rcpBottomRow, '.png'), width=width, height=height, res=450)
			
				# # # pars <- par(mfrow=c(4, 4), oma=rep(0, 4), mar=c(0, 0, 0, 0))

				# # # # by EACH RASTER
				# # # for (map in maps) {
					
					# # # say('subplot ', map)
					
					# # # lab <- if (is.na(map)) {
						# # # NA
					# # # } else if (map == 'topRow') {
						# # # 'Current climate'
					# # # } else if (map == 'middleRow') {
						# # # paste0(futYear, ' climate RCP ', rcpMiddleRowNice)
					# # # } else if (map == 'bottomRow') {
						# # # paste0(futYear, ' climate RCP ', rcpBottomRowNice)
					# # # } else if (map == 'leftCol') {
						# # # '2014 forest'
					# # # } else if (map == 'middleCol') {
						# # # paste0(futYear, ' forest\nstrict', ifelse(inset == 3, '\n', ' '), 'protection')
					# # # } else if (map == 'rightCol') {
						# # # paste0(futYear, ' forest\nrelaxed', ifelse(inset == 3, '\n', ' '), 'protection')
					# # # }
					
					# # # rowCol <- if (is.na(map)) { c(1, 1) } else
						# # # if (map == 'topRow') { c(2, 1)} else
						# # # if (map == 'middleRow') { c(3, 1)} else
						# # # if (map == 'bottomRow') { c(4, 1)} else
						# # # if (map == 'leftCol') { c(1, 2)} else
						# # # if (map == 'middleCol') { c(1, 3)} else
						# # # if (map == 'rightCol') { c(1, 4)} else
						# # # if (map == 'r1c1') { c(1, 1) + 1 } else
						# # # if (map == 'r2c1') { c(1, 2) + 1 } else
						# # # if (map == 'r3c1') { c(1, 3) + 1 } else
						# # # if (map == 'r1c2') { c(2, 1) + 1 } else
						# # # if (map == 'r2c2') { c(2, 2) + 1 } else
						# # # if (map == 'r3c2') { c(2, 3) + 1 } else
						# # # if (map == 'll') { c(3, 1) + 1 } else
						# # # if (map == 'lm') { c(3, 2) + 1 } else
						# # # if (map == 'lr') { c(3, 3) + 1 }
					
					# # # letter <- if (is.na(map)) { NA } else
						# # # if (map == 'r1c1') { 'a' } else
						# # # if (map == 'r2c1') { 'b' } else
						# # # if (map == 'r3c1') { 'c' } else
						# # # if (map == 'r1c2') { 'd' } else
						# # # if (map == 'r2c2') { 'e' } else
						# # # if (map == 'r3c2') { 'f' } else
						# # # if (map == 'll') { 'g' } else
						# # # if (map == 'lm') { 'h' } else
						# # # if (map == 'lr') { 'i' }
					
					# # # par(mfg=rowCol)
					# # # plot.new()
					
					# # # # top left corner
					# # # if (is.na(map)) {
					
						# # # 'hey!'
					
					# # # # column labels
					# # # } else if (map %in% c('leftCol', 'middleCol', 'rightCol')) {

						# # # plot.window(xlim=c(-1, 1), ylim=c(-1, 1))
						# # # if (inset == 1 | inset == 2) text(0, -0.70, labels=lab, cex=1.4, xpd=NA, pos=1, font=2)
						# # # if (inset == 3) text(0, -0.50, labels=lab, cex=1.4, xpd=NA, pos=1, font=2)
						
					# # # # row labels
					# # # } else if (map %in% c('topRow', 'middleRow', 'bottomRow')) {

						# # # plot.window(xlim=c(-1, 1), ylim=c(-1, 1))
						# # # text(0.90, 0, labels=lab, cex=1.4, xpd=NA, srt=90, font=2)
						
					# # # } else {
						
						# # # x <- get(map)
						# # # xCrop <- crop(x, focus)
						# # # vals <- cellStats(xCrop, 'sum')
						# # # insetCols <- if (vals == 0) { cols[1] } else { cols }

						# # # # focal plot
						# # # plot(focus, ann=FALSE, xaxt='s')
						# # # # plot(hsCrop, col=hsCols, add=TRUE, legend=FALSE)
						# # # # plot(humidForestBufferMaskCrop_utm38s, col='gray80', legend=FALSE, bty='n', xaxt='n', yaxt='n', fg='white', add=TRUE)
						# # # if (!is.null(pasCrop)) plot(pasCrop, border='blue', col=alpha('blue', 0.3), add=TRUE)
						# # # plot(xCrop, col=cols, breaks=seq(0, 100, length.out=colBreaks), legend=FALSE, add=TRUE)
						# # # if (!is.null(pasCrop)) plot(pasCrop, border='blue', col=alpha('blue', 0.1), add=TRUE)
						# # # plot(madagascarCrop_utm38s, border='black', xpd=NA, add=TRUE)

						# # # xRange <- ext@xmax - ext@xmin
						# # # yRange <- ext@ymax - ext@ymin

						# # # # scale bar
						# # # if (map == 'lr') {
							
							# # # scale <- 25000 # meters
							
							# # # xStart <- ext@xmax - 1.25 * scale
							# # # yAt <- 0.03 * yRange + ext@ymin
							
							# # # lines(c(xStart, xStart + scale), c(yAt, yAt), lwd=4, lend=2, xpd=NA)
							# # # text(xStart + 0.5 * scale, yAt + 0.05 * yRange, labels=paste(scale / 1000, 'km'), cex=1.2)
							
						# # # }

						# # # # title
						# # # usr <- par('usr')
						# # # titleX <- 0.01 * xRange + ext@xmin
						# # # titleY <- 0.92 * yRange + ext@ymin
						
						# # # text(titleX, titleY, labels=paste0(letter, ')'), pos=4, cex=1.3, xpd=NA, font=2)

					# # # } # if map
						
				# # # } # next map
				
				# # # title(main=date(), outer=TRUE, cex.sub=1, line=-2)
				
			# # # dev.off()
			
		# # # } # next inset
				
		# # # par(pars)
	
	
# say('################################################################')
# say('### create 3D displays of ecological niche model predictions ###')
# say('################################################################')

	# ### definitions
	# ###############

	# # futYear <- '2050'
	# futYear <- '2070'

	# z <- 7.5

	# # ENM rasters
	# upperLeft <- raster(paste0('./Ecological Niche Models/Prediction Rasters/glmEnm_vareciaGenus_climateCurrent_forest2014_utm38s.tif'))
	# upperMiddle <- raster(paste0('./Ecological Niche Models/Prediction Rasters/glmEnm_vareciaGenus_climateCurrent_forest', futYear, '_defoOutsidePAs_utm38s.tif'))
	# upperRight <- raster(paste0('./Ecological Niche Models/Prediction Rasters/glmEnm_vareciaGenus_climateCurrent_forest', futYear, '_defoAnywhere_utm38s.tif'))
	
	# middleLeft <- raster(paste0('./Ecological Niche Models/Prediction Rasters/glmEnm_vareciaGenus_climate2070EnsembleMeanRcp4pt5_forest2014_utm38s.tif'))
	# middleMiddle <- raster(paste0('./Ecological Niche Models/Prediction Rasters/glmEnm_vareciaGenus_climate', futYear, 'EnsembleMeanRcp4pt5_forest', futYear, '_defoOutsidePAs_utm38s.tif'))
	# middleRight <- raster(paste0('./Ecological Niche Models/Prediction Rasters/glmEnm_vareciaGenus_climate', futYear, 'EnsembleMeanRcp4pt5_forest', futYear, '_defoAnywhere_utm38s.tif'))
	
	# lowerLeft <- raster(paste0('./Ecological Niche Models/Prediction Rasters/glmEnm_vareciaGenus_climate', futYear, 'EnsembleMeanRcp8pt5_forest2014_utm38s.tif'))
	# lowerMiddle <- raster(paste0('./Ecological Niche Models/Prediction Rasters/glmEnm_vareciaGenus_climate', futYear, 'EnsembleMeanRcp8pt5_forest', futYear, '_defoOutsidePAs_utm38s.tif'))
	# lowerRight <- raster(paste0('./Ecological Niche Models/Prediction Rasters/glmEnm_vareciaGenus_climate', futYear, 'EnsembleMeanRcp8pt5_forest', futYear, '_defoAnywhere_utm38s.tif'))

	# panels <- c('upperLeft', 'upperMiddle', 'upperRight', 'middleLeft', 'middleMiddle', 'middleRight', 'lowerLeft', 'lowerMiddle', 'lowerRight')

	# # elevation
	# elev <- raster('./Data/Topography - GMTED2010/elevationGmted2010_utm38s.tif')
	# projection(elev) <- madEaProj
	
	# # # occurrences
	# # load('./Ecological Niche Models/Collated Presence and Background Point Data for Varecia VARIEGATA.RData')
	# # vv <- taxonData[taxonData$presBg == 1, ]
	# # load('./Ecological Niche Models/Collated Presence and Background Point Data for Varecia RUBRA.RData')
	# # vr <- taxonData[taxonData$presBg == 1, ]
	
	# # vv <- SpatialPoints(vv[ , longLat], getCRS('wgs84', TRUE))
	# # vr <- SpatialPoints(vr[ , longLat], getCRS('wgs84', TRUE))
	
	# # vv <- sp::spTransform(vv, CRS(madEaProj))
	# # vr <- sp::spTransform(vr, CRS(madEaProj))
	
	# # vv <- elimCellDups(vv, elev)
	# # vr <- elimCellDups(vr, elev)
	
	# # vv <- buffer(vv, width=360)
	# # vr <- buffer(vr, width=360)
	
	# ### by inset
	# ############
	
	# for (inset in 1:3) {

		# say(inset, level=2)
	
		# focus <- get(paste0('madFocus', inset))
		# elevCrop <- crop(elev, focus)
		
		# # ### occurrences
		# # thisVv <- crop(vv, elev)
		# # thisVr <- crop(vr, elev)
		
		# # if (!is.null(thisVv)) {
			# # vvRast <- rasterize(thisVv, elev, 1, background=0)
			# # vvMat <- as.matrix(vvRast)
		# # }

		# # if (!is.null(thisVr)) {
			# # vrRast <- rasterize(thisVr, elev, 1, background=0)
			# # vrMat <- as.matrix(vrRast)
		# # }
		
		# ### hillshading
		# elevMat <- matrix(
			# raster::extract(elevCrop, raster::extent(elevCrop), method = 'bilinear'),
			# nrow = ncol(elevCrop),
			# ncol = nrow(elevCrop)
		# )

		# hs <- elevMat %>%
			# sphere_shade(texture = 'imhof4') %>%
			# add_shadow(ambient_shade(elevMat, zscale=z, anglebreaks = seq(65, 65, 1))) %>%
			# add_shadow(ray_shade(elevMat, zscale=z, lambert=FALSE, anglebreaks = seq(65, 65, 1)))
		
		# ### protected areas
		# # define boundaries as difference between inner and outer buffer
		# load('./Data/Protected Areas/WDPA_Sept2018_MDG-shapefile-polygons-onlyTerrestrial.RData')
		# pas <- crop(pas, focus)

		# inner <- gBuffer(pas, width=-150)
		# outer <- gBuffer(pas, width=150)

		# inner <- rasterize(inner, elevCrop, 1)
		# outer <- rasterize(outer, elevCrop, 1)

		# naTo0 <- function(x) ifelse(is.na(x), 0, x)

		# beginCluster(4)
			# inner <- clusterR(inner, calc, args=list(fun=naTo0))
			# outer <- clusterR(outer, calc, args=list(fun=naTo0))
		# endCluster()

		# paBorders <- outer - inner
		
		# paBorders <- matrix(
			# raster::extract(paBorders, raster::extent(paBorders)),
			# nrow = ncol(paBorders),
			# ncol = nrow(paBorders)
		# )
		
		# paBorders <- t(paBorders)
		
		# if (anyNA(paBorders)) paBorders[is.na(paBorders)] <- 0

		# ### by enm
		# ##########
		
		# for (panel in panels) {
		
			# say('inset ', inset, ' ', panel)
		
			# enm <- get(panel)
			# enm <- crop(enm, focus)
			# enm <- enm / 100
			
			# # convert enm to 3 layers RGB
			# r <- matrix(
				# raster::extract(enm, raster::extent(enm), method = 'bilinear'),
				# nrow = ncol(enm),
				# ncol = nrow(enm)
			# )
			# r <- enm
			# g <- 0 * r
			# b <- 0 * r
			# alpha <- r

			# r <- as.matrix(r)
			# g <- as.matrix(g)
			# b <- as.matrix(b)
			# alpha <- as.matrix(alpha)

			# if (anyNA(r)) r[is.na(r)] <- 0
			# if (anyNA(g)) g[is.na(g)] <- 0
			# if (anyNA(b)) b[is.na(b)] <- 0
			# if (anyNA(alpha)) alpha[is.na(alpha)] <- 0

			# alpha <- 0.8 * alpha

			# rgb <- array(c(r, g, b, alpha), c(nrow(elevCrop), ncol(elevCrop), 4))

			# hsEnm <- add_overlay(hs, rgb, 4)

			# ## add PA boundaries
			# hsEnmPas <- hsEnm
			# for (i in 1:nrow(paBorders)) {
				# for (j in 1:ncol(paBorders)) {
					# if (paBorders[i, j] == 1) {
						# hsEnmPas[i, j, 1] <- 0
						# hsEnmPas[i, j, 2] <- 0
						# hsEnmPas[i, j, 3] <- 1
					# }
				# }
			# }

			# # ### add occurrences
			
			# # # variegata
			# # if (!is.null(thisVv)) {
			
				# # for (i in 1:nrow(vvMat)) {
					# # for (j in 1:ncol(vvMat)) {
						# # if (vvMat[i, j] == 1) {
							# # # hsEnmPas[i, j, 1] <- 0
							# # # hsEnmPas[i, j, 2] <- 248 / 255
							# # # hsEnmPas[i, j, 3] <- 252 / 255
							# # hsEnmPas[i, j, 1] <- 0
							# # hsEnmPas[i, j, 2] <- 0
							# # hsEnmPas[i, j, 3] <- 0
						# # }
					# # }
				# # }
				
			# # }
				
			# # # rubra
			# # if (!is.null(thisVr)) {
			
				# # for (i in 1:nrow(vrMat)) {
					# # for (j in 1:ncol(vrMat)) {
						# # if (vrMat[i, j] == 1) {
							# # # hsEnmPas[i, j, 1] <- 252 / 255
							# # # hsEnmPas[i, j, 2] <- 248 / 255
							# # # hsEnmPas[i, j, 3] <- 0
							# # hsEnmPas[i, j, 1] <- 0
							# # hsEnmPas[i, j, 2] <- 0
							# # hsEnmPas[i, j, 3] <- 0
						# # }
					# # }
				# # }
				
			# # }
			
			# ### plot!
			# solidDepth <- min(elevMat, na.rm=TRUE)
			# plot_3d(hsEnmPas, elevMat, zscale=z, theta=45, phi=30, water=TRUE, zoom=0.8, fov=60, soliddepth=solidDepth, windowsize=c(1600, 1200), shadow=FALSE, background='white')
			
			# outDir <- paste0('./Figures & Tables/Ecological Niche Models - Predictions 3D/Panels for ', futYear)
			# dirCreate(outDir)
			
			# rgl.snapshot(paste0(outDir, '/Inset ', inset, ' Panel ', panel, '.png'), fmt='png', top=FALSE)
			# rgl.close()

		# } # next panel/ENM
		
	# } # next inset

	# ### create scale bar
	# ####################
	
	# png('./Figures & Tables/Ecological Niche Models - Predictions 3D/Legend for ENM.png', width=800, height=800, res=300)
		
		# plot(0, 0, col='white', fg='white', ann=FALSE, xaxt='n', yaxt='n')
		# paSwatch <- list(swatchAdjY=c(0, 0.12), col=NA, border=alpha('blue', 1), labels='PA')
		# legendGrad('left', inset=0.4, height=1.8, labels=c(0, 0.5, 1), labAdjX=0.6, col=c('white', 'white', alpha('red', 0.2), 'red'), title='Suitability', adjX=c(0.1, 0.7), adjY=c(0.25, 0.75), boxBorder=NA, swatches=list(paSwatch), cex=1, xpd=NA)
		
	# dev.off()
		
# say('############################################')
# say('### create 3D displays of climate change ###')
# say('############################################')

	# say('Having problems with colorizing climate... probably related to resampling of climate rasters from 30 arcsec to 30 m resolution.')

	# ### definitions
	# ###############

	# # futYear <- '2050'
	# futYear <- '2070'

	# z <- 800 / 4

	# humidForestBufferMask_utm38s <- raster('./Study Region & Masks/UTM 38S 30-m Resolution/humidForestBufferMask_utm38s.tif')
	# load('./Study Region & Masks/WGS84 30-arcsec Resolution/Eastern Humid Forest Polygon Buffer.RData')
	
	# # climate rasters
	# current <- raster('./Data/Climate - WORLDCLIM Ver 2 Rel 1 Jun 2016 (1970-2000)/bio10.tif')
	# current <- crop(current, humidForestBuffer_wgs84)

	# panels <- c('current')

	# # elevation
	# elev <- raster('./Data/Topography - WORLDCLIM Ver 1pt4 Rel 3/elevation_wgs84.tif')
	
	# ### by inset
	# ############
	
	# # for (inset in 1:3) {
	# for (inset in 3) {

		# say(inset, level=2)
	
		# focus_utm38s <- get(paste0('madFocus', inset))
		# focus_wgs84 <- sp::spTransform(focus_utm38s, getCRS('wgs84', TRUE))
		# elevCrop <- crop(elev, focus_wgs84)
		
		# ### hillshading
		# elevMat <- matrix(
			# raster::extract(elevCrop, raster::extent(elevCrop), method = 'bilinear'),
			# nrow = ncol(elevCrop),
			# ncol = nrow(elevCrop)
		# )

		# hs <- elevMat %>%
			# sphere_shade(texture = 'imhof4') %>%
			# add_shadow(ambient_shade(elevMat, zscale=z, anglebreaks = seq(65, 65, 1))) %>%
			# add_shadow(ray_shade(elevMat, zscale=z, lambert=FALSE, anglebreaks = seq(65, 65, 1)))

		# ### protected areas
		# # define boundaries as difference between inner and outer buffer
		# load('./Data/Protected Areas/WDPA_Sept2018_MDG-shapefile-polygons-onlyTerrestrial.RData')
		# pas <- crop(pas, focus_utm38s)

		# inner <- gBuffer(pas, width=-400)
		# outer <- gBuffer(pas, width=400)

		# inner <- sp::spTransform(inner, getCRS('wgs84', TRUE))
		# outer <- sp::spTransform(outer, getCRS('wgs84', TRUE))
		
		# inner <- rasterize(inner, elevCrop, 1)
		# outer <- rasterize(outer, elevCrop, 1)

		# naTo0 <- function(x) ifelse(is.na(x), 0, x)

		# inner <- calc(inner, fun=naTo0)
		# outer <- calc(outer, fun=naTo0)

		# paBorders <- outer - inner
		
		# paBorders <- matrix(
			# raster::extract(paBorders, raster::extent(paBorders)),
			# nrow = ncol(paBorders),
			# ncol = nrow(paBorders)
		# )
		
		# paBorders <- t(paBorders)
		
		# if (anyNA(paBorders)) paBorders[is.na(paBorders)] <- 0

		# ### by climate scenario
		# #######################
		
		# # NB "x" is the "color" raster (climate in this case)

		# focusBuffUtm38s <- buffer(focus_utm38s, 5000)
		# focusBuff_wgs84 <- sp::spTransform(focusBuffUtm38s, getCRS('wgs84', TRUE))
		
		# for (panel in panels) {
		
			# say('inset ', inset, ' ', panel)
		
			# x <- get(panel)
			# elevCrop <- crop(elev, focus_wgs84)
			
			# ### convert x to 3 layers RGB
			# clim <- matrix(
				# raster::extract(x, raster::extent(x), method = 'bilinear'),
				# nrow = ncol(x),
				# ncol = nrow(x)
			# )
			
			# r <- g <- b <- NA * clim
			
			# # # match palette colors to climate values
			# # # pal <- c('blue', 'yellow', 'red')
			# # # pal <- c('darkblue', 'blue', 'cyan', 'yellow', 'orange', 'red', 'darkred')
			# # # pal <- c('blue', 'red')
			# # pal <- c('blue', 'yellow')
			# # pal <- colorRampPalette(pal)
			# # cols <- pal(40)
			
			# # rgbCols <- col2rgb(cols)
			# # rgbCols <- rgbCols / 255
			
			# # minClim <- minValue(x)
			# # maxClim <- maxValue(x)
			
			# # cols <- data.frame(col=cols, r=rgbCols[1, ], g=rgbCols[2, ], b=rgbCols[3, ], climate=seq(minClim, maxClim, length.out=length(cols)))
			
			# # for (row in 1:nrow(clim)) {
				# # for (col in 1:ncol(clim)) {
					# # index <- which.min(abs(cols$climate - clim[row, col]))
					# # if (length(index) > 0) {
						# # r[row, col] <- cols$r[index]
						# # g[row, col] <- cols$g[index]
						# # b[row, col] <- cols$b[index]
					# # }
				# # }
			# # }

			# r <- 1 + (0 * clim)
			# g <- b <- 0 * r
			
			# # alpha <- 0.6 + (0 * clim)
			# alpha <- 1 + (0 * clim)
			# alpha <- as.matrix(x - minValue(x))
			# alpha <- alpha / max(alpha)
			# alpha <- t(alpha)
			
			# if (anyNA(r)) r[is.na(r)] <- 0
			# if (anyNA(g)) g[is.na(g)] <- 0
			# if (anyNA(b)) b[is.na(b)] <- 0
			# if (anyNA(alpha)) alpha[is.na(alpha)] <- 0

			# rgb <- array(c(r, g, b, alpha), c(nrow(elevCrop), ncol(elevCrop), 4))

			# hsClim <- add_overlay(hs, rgb, 4)
			
# solidDepth <- min(elevMat, na.rm=TRUE)
# plot_3d(hsClim, elevMat, zscale=z, theta=45, phi=30, water=TRUE, zoom=0.8, fov=60, soliddepth=solidDepth, windowsize=c(1600, 1200), shadow=FALSE, background='white')			
			
			# ## add PA boundaries
			# hsClimPas <- hsClim
			# for (i in 1:nrow(paBorders)) {
				# for (j in 1:ncol(paBorders)) {
					# if (paBorders[i, j] == 1) {
						# hsClimPas[i, j, 1] <- 0
						# hsClimPas[i, j, 2] <- 0
						# hsClimPas[i, j, 3] <- 1
					# }
				# }
			# }

			
			# ### plot!
			# solidDepth <- min(elevMat, na.rm=TRUE)
			# plot_3d(hsClimPas, elevMat, zscale=z, theta=45, phi=30, water=TRUE, zoom=0.8, fov=60, soliddepth=solidDepth, windowsize=c(1600, 1200), shadow=FALSE, background='white')
			
			# outDir <- paste0('./Figures & Tables/Climate - 3D Maps/Panels for ', futYear)
			# dirCreate(outDir)
			
			# rgl.snapshot(paste0(outDir, '/Inset ', inset, ' Panel ', panel, '.png'), fmt='png', top=FALSE)
			# rgl.close()

		# } # next panel/ENM
		
	# } # next inset

	# ### create scale bar
	# ####################
	
	# # # png('./Figures & Tables/Ecological Niche Models - Predictions 3D/Legend for ENM.png', width=800, height=800, res=300)
		
		# # # plot(0, 0, col='white', fg='white', ann=FALSE, xaxt='n', yaxt='n')
		# # # paSwatch <- list(swatchAdjY=c(0, 0.12), col=NA, border=alpha('blue', 1), labels='PA')
		# # # legendGrad('left', inset=0.4, height=1.8, labels=c(0, 0.5, 1), labAdjX=0.6, col=c('white', 'white', alpha('red', 0.2), 'red'), title='Suitability', adjX=c(0.1, 0.7), adjY=c(0.25, 0.75), boxBorder=NA, swatches=list(paSwatch), cex=1, xpd=NA)
		
	# # # dev.off()
		
# say('#############################################')
# say('### create display range maps by surveyor ###')
# say('#############################################')

	# # generalization
	# outDir <- './Figures & Tables/Maps of Occurrences/'
	# dirCreate(outDir)

	# # ancillary geo data
	# elev <- raster('./Data/Topography - GMTED2010/elevationGmted2010_utm38s.tif')
	
	# load('./Study Region & Masks/UTM 38S 30-m Resolution/Eastern Humid Forest Polygon.RData')
	# load('./Study Region & Masks/UTM 38S 30-m Resolution/Madagascar from GADM 3.6.RData')

	# # taxon data
	# varecia <- read.csv('./Data/Varecia/05 Varecia Occurrences with Weights Based on Neighbor Proximity.csv')
	# varecia <- SpatialPointsDataFrame(varecia[ , longLat], data=varecia, proj4=getCRS('wgs84', TRUE))
	# varecia <- sp::spTransform(varecia, CRS(madEaProj))
	
	# # # elevation background
	# # cols <- colorRampPalette(c('black', 'gray50', 'white'))
	# # cols <- cols(20)
	# # for (i in seq_along(cols)) cols[i] <- alpha(cols[i], 0.5)
	
	# ### Madagascar
	# ##############
		
		# providers <- sort(unique(varecia$provider))
			
		# for (i in seq_along(providers)) {
		
			# provider <- providers[i]
			
			# providerNice <- gsub(provider, pattern='_', replacement=' & ')
			# if (providerNice == 'BrownAndYoder') providerNice <- 'Brown & Yoder'
		
			# png(paste0(outDir, '/Varecia Surveys by ', providerNice, '.png'), width=1000, height=1400, res=450)
		
				# pars <- par(oma=c(1, 1, 1, 1), mar=rep(0.1, 4))

				# # main plot
				# plot(humidForest_utm38s, bty='n', ann=FALSE, col=alpha('lightgreen', 0.4), main=NA)
				# plot(madagascar_utm38s, lwd=0.2, add=TRUE)

				# # survey sites
				# vAbs <- varecia[varecia$provider == provider & varecia$presAbs == 0, ]
				# vvPres <- varecia[varecia$species == 'Varecia variegata' & varecia$provider == provider & varecia$presAbs == 1, ]
				# vrPres <- varecia[varecia$species == 'Varecia rubra' & varecia$provider == provider & varecia$presAbs == 1, ]
				
				# # remove identical points among pres/abs (presence sites for one species were often included as absences for the other)
				# dups <- integer()
				# absCoords <- coordinates(vAbs)
				# vrPresCoords <- coordinates(vrPres)
				# vvPresCoords <- coordinates(vvPres)
				
				# if (nrow(absCoords) > 0) {
					
					# for (i in 1:nrow(absCoords)) {
					
						# if (nrow(vvPresCoords) > 0) {
					
							# if (absCoords[i, 1] %in% vvPresCoords[ , 1] & absCoords[i, 2] %in% vvPresCoords[ , 2]) dups <- c(dups, i)							
						# }
						
						# if (nrow(vrPresCoords) > 0) {
					
							# if (absCoords[i, 1] %in% vrPresCoords[ , 1] & absCoords[i, 2] %in% vrPresCoords[ , 2]) dups <- c(dups, i)
							
						# }
						
					# }
					
					# if (length(dups) > 0) vAbs <- vAbs[-dups, ]
					
				# }
				
				# if (length(vAbs) > 0) points(vAbs, pch=1, col='red', cex=0.4)
				# if (length(vvPres) > 0) points(vvPres, pch=3, col='black', cex=0.4)
				# if (length(vrPres) > 0) points(vrPres, pch=4, col='black', cex=0.4)
			
				# legend(
					# 'bottomright',
					# inset=c(0, 0.15),
					# bty='n',
					# legend=c('V. variegata detection', 'V. rubra detection', 'Non-detection'),
					# col=c('black', 'black', 'red'),
					# pt.bg=c('chartreuse', 'chartreuse', NA),
					# pch=c(3, 4, 1),
					# cex=0.4
				# )
					
				# title(main=providerNice, cex.main=0.8, line=-0.5)
				# # title(sub=date(), outer=TRUE, cex.sub=0.3, line=0)
				
			# dev.off()

		# }
			
		# par(pars)
		
# say('##############################')
# say('### niche overlap analysis ###')
# say('##############################')

	# load('./Ecological Niche Models/Collated Presence and Background Point Data for Varecia RUBRA.RData')
	# vr <- taxonData[taxonData$presBg == 1, longLat]
	# load('./Ecological Niche Models/Collated Presence and Background Point Data for Varecia VARIEGATA.RData')
	# vv <- taxonData[taxonData$presBg == 1, longLat]
	
	# mask_wgs84 <- raster('./Study Region & Masks/WGS84 30-arcsec Resolution/mask_wgs84.tif')

	# vv <- SpatialPoints(vv, getCRS('wgs84', TRUE))
	# vr <- SpatialPoints(vr, getCRS('wgs84', TRUE))
	
	# vr <- elimCellDups(vr, mask_wgs84)
	# vv <- elimCellDups(vv, mask_wgs84)
	
	# vr <- geoThin(vr, minDist=70000 / 8)
	# vv <- geoThin(vv, minDist=70000 / 4)
	
	# tol1 <- tol2 <- tol12 <- 70000 / 8
	
	# rands <- randPointsBatch(
		# 'randPointsRespectingSelfOther2',
		# iterations=100,
		# x1=vv,
		# x2=vr,
		# tol1=tol1,
		# tol2=tol2,
		# tol12=tol12,
		# rast=mask_wgs84,
		# verbose=TRUE,
		# verboseEach=FALSE
	# )

	# # match climate to sampled sites
	# clim <- raster::stack(paste0('./Data/Climate - WORLDCLIM Ver 2 Rel 1 Jun 2016 (1970-2000)/bio', prefix(bioclims, 2), '.tif'))
	# rands <- randPointsBatchExtract(rands, clim, TRUE)

	# # remove NAs
	# for (i in 1:100) {
		
		# x1 <- rands[[i]]$x1rand
		# x2 <- rands[[i]]$x2rand
		
		# nas1 <- naRows(x1@data)
		# nas2 <- naRows(x2@data)
		
		# if (length(nas1) > 0) x1 <- x1[-nas1, ]
		# if (length(nas2) > 0) x2 <- x2[-nas2, ]
		
		# rands[[i]]$x1rand <- x1
		# rands[[i]]$x2rand <- x2
		
	# }
	
	# sampled <- randPointsBatchSampled(rands)
	
	# # create common PCA of environmental space
	# bgPca <- princomp(sampled@data[ , paste0('bio', prefix(bioclims, 2))], cor=TRUE)

	# vrEnv <- extract(clim, vr)
	# vvEnv <- extract(clim, vv)

	# # observed niche overlap
	# obsOverlap <- enmSdm::nicheOverlap(
		# x1=vrEnv,
		# x2=vvEnv,
		# env=bgPca,
		# vars=paste0('bio', prefix(bioclims, 2)),
		# bins=100,
		# cor=TRUE
	# )
	
	# # null niche overlap
	# nullOverlap <- randPointsBatchNicheOverlap(
		# rands=rands,
		# env=bgPca,
		# vars=paste0('bio', prefix(bioclims, 2)),
		# bins=100,
		# cor=TRUE
	# )
	
	# outDir <- './Niche Overlap Analysis'
	# dirCreate(outDir)
	# save(rands, file=paste0(outDir, '/Randomized Point Locations with Environment.RData'))
	# write.csv(obsOverlap, paste0(outDir, '/Observed Niche Overlap.csv'), row.names=FALSE)
	# write.csv(nullOverlap, paste0(outDir, '/Null Niche Overlap.csv'), row.names=FALSE)

	# # report niche overlap statistics
	# dObs <- obsOverlap[['d']]
	# rankCorObs <- obsOverlap[['rankCor']]
	
	# dQuant <- sum(nullOverlap$d <= dObs) / 100
	# rankCorQuant <- sum(nullOverlap$rankCor <= rankCorObs) / 100
	
	# sink('./Niche Overlap Analysis/Niche Overlap Analysis Results.txt', split=TRUE)
		
		# say('### niche overlap analysis between Varecia variegata and V. rubra ###')
		# say('')
		# say('=== Warren\'s D: ===')
		# say('Observed D:', dObs)
		# say('P value:', 1 - dQuant)
		# say('')
		# say('=== rank correlation: ===')
		# say('Observed Rank Correlation: ', rankCorObs)
		# say('P value:', 1 - rankCorQuant)

	# sink()	
	
# # # say('#############################################################')
# # # say('### analysis of spatial autocorrelation within predictors ###')
# # # say('#############################################################')

	# # # say('I am calculating the characteristic scale of spatial autocorrelation within each predictor as the distance where the median of the observed distribution of environmental distances falls at or above the median of the distribution of environmental distances using scrambled environmental distances. Samples are from points randomly located across Madagascar.', breaks=80)
	
	# # # ### measure characteristic scale of spatial autocorrelation across predictors
	# # # #############################################################################

		# # # bioPredictors <- c(paste0('bio', prefix(bioclims, 2)))
		
		# # # load('./Ecological Niche Models/Random Background Sites from across Madagascar with Current Environmental Data.RData')
		
		# # # # calculate geographic distance between sites
		# # # geoDists <- pointDist(bg[ , longLat])
		# # # diag(geoDists) <- NA
		# # # geoDists <- geoDists / 1000
		# # # geoDists <- c(geoDists)
		
		# # # # stores characteristic spatial distances for bioPredictors
		# # # charDists <- data.frame()
		
		# # # for (thisPred in bioPredictors) {

			# # # say(thisPred)
		
			# # # # calculate environmental distances between sites
			# # # envDists <- matrix(NA, nrow=nrow(bg), ncol=nrow(bg))
			
			# # # for (i in 1:nrow(bg)) envDists[i, ] <- abs(bg[i, thisPred] - bg[ , thisPred])
			# # # envDists <- c(envDists)
			
			# # # # bin geographic distances
			# # # geoDistStepSize <- 20
			# # # maxGeoDist <- max(geoDists, na.rm=TRUE)
			# # # maxGeoDistNice <- geoDistStepSize * ceiling(maxGeoDist / geoDistStepSize)
			
			# # # geoBins <- seq(0, maxGeoDistNice, by=geoDistStepSize)
			# # # envFreqLow <- envFreqMedian <- envFreqHigh <- rep(0, length(geoBins))
			# # # envFreqLowRand <- envFreqMedianRand <- envFreqHighRand <- rep(0, length(geoBins))

			# # # # tally distribution of environmental distances in each distance bin
			# # # geoFreq <- rep(0, length(geoBins))
			# # # for (i in 1:(length(geoBins) - 1)) {
				
				# # # # tally geographic distance in bin
				# # # whichGeoDistsInBin <- which(geoDists >= geoBins[i] & geoDists < geoBins[i + 1])
				# # # geoFreq[i] <- length(whichGeoDistsInBin)
				
				# # # # tally environmental distance in bin
				# # # envDistsIn <- envDists[whichGeoDistsInBin]
				
				# # # if (length(envDistsIn) > 0) {
					
					# # # envQuants <- quantile(envDistsIn, c(0.025, 0.5, 0.975))
					
					# # # envFreqLow[i] <- envQuants[1]
					# # # envFreqMedian[i] <- envQuants[2]
					# # # envFreqHigh[i] <- envQuants[3]
					
					# # # envDistsInRand <- sample(envDists, length(whichGeoDistsInBin))

					# # # envQuantsRand <- quantile(envDistsInRand, c(0.025, 0.5, 0.975))
					
					# # # envFreqLowRand[i] <- envQuantsRand[1]
					# # # envFreqMedianRand[i] <- envQuantsRand[2]
					# # # envFreqHighRand[i] <- envQuantsRand[3]
						
				# # # }
					
			# # # }
			
			# # # # plot
			# # # dirCreate('./Analysis - Global Change Vulnerability')
			
			# # # png(paste0('./Analysis - Global Change Vulnerability/Characteristic Scale of Spatial Autocorrelation for ', thisPred, '.png'), width=1200, height=600, res=150)
				
				# # # yMax <- 10 * ceiling(max(envDists / 10, na.rm=TRUE))
				# # # plot(geoBins + geoDistStepSize, geoBins, ylim=c(0, yMax), col='white', xlab='Distance (km)', ylab='Environmental distance', main=thisPred)
				
				# # # nudge <- 0.01 * geoDistStepSize
				
				# # # for (i in seq_along(geoBins)) {
					
					# # # if (geoFreq[i] > 0) {
					
						# # # x <- geoBins[i] + 0.5 * geoDistStepSize
						# # # lines(c(x, x), c(envFreqLow[i], envFreqHigh[i]), lwd=8, lend=2, col='gray')
						# # # points(x, envFreqMedian[i], pch=15)
						# # # lines(c(x, x) + nudge, c(envFreqLowRand[i], envFreqHighRand[i]), lwd=1, col='red')
						# # # points(x + nudge, envFreqMedianRand[i], pch=15, col='darkred')
					
					# # # }
				# # # }
				
				# # # # at what geographic distance does median of observed environmental distances surpass randomized environmental distance?
				# # # whichCharGeoDist <- which(envFreqMedian >= envFreqMedianRand)[1]
				# # # charGeoDist <- geoBins[whichCharGeoDist]

				# # # y0 <- 0.12 * yMax + envFreqHighRand[whichCharGeoDist]
				# # # y1 <- 0.01 * yMax + envFreqHighRand[whichCharGeoDist]
				
				# # # arrows(x0=charGeoDist + 0.5 * geoDistStepSize, x1=charGeoDist + 0.5 * geoDistStepSize, y0=y0, y1=y1, angle=17, length=0.1)

				# # # legend('topleft', inset=0.01, legend=c('Observed', 'Random'), col=c('gray', 'red'), pch=c(15, 15), lwd=c(4, 1), bty='n')
				
				# # # # remember
				# # # thisCharDist <- data.frame(
					# # # pred=thisPred,
					# # # charDist_km=charGeoDist - 0.5 * geoDistStepSize
				# # # )
				
				# # # charDists <- rbind(charDists, thisCharDist)
				
			# # # dev.off()
			
		# # # } # next predictor
	
		# # # write.csv(charDists, './Analysis - Global Change Vulnerability/Characteristic Scale of Spatial Autocorrelation for each ENM Predictor.csv', row.names=FALSE)
	
say('DONE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ', date(), level=1, deco='%')
