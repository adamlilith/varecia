### Modeling range of Varecia, the Black-and-White Lemur (with Toni Lyn Morelli)
### Adam B. Smith | Missouri Botanical Garden | adamDOTsmithATmobotDOTorg | 2018-02

### This code is version 6 of attempts to model the exposure of two lemurs, Varecia variegata and Varecia rubra, to anticipated climate change and deforestation in Madagascar. The analysis depends on two models, one of deforestation and an ecological niche model using forest cover and climate as predictors.

# source('C:/Ecology/Drive/Research/Varecia (Toni Lyn Morelli)/Versions 06/Code/Varecia Ver 06 Forest Modeling.r')
# source('E:/Ecology/Drive/Research/Varecia (Toni Lyn Morelli)/Versions 06/Code/Varecia Ver 06 Forest Modeling.r')
# source('H:/Global Change Program/Research/Varecia (Toni Lyn Morelli)/Versions 06/Code/Varecia Ver 06 Forest Modeling.r')

	setwd('C:/ecology/Drive/Research/Varecia (Toni Lyn Morelli)/Versions 06')
	# setwd('E:/ecology/Drive/Research/Varecia (Toni Lyn Morelli)/Versions 06')
	# setwd('H:/Global Change Program/Research/Varecia (Toni Lyn Morelli)/Versions 06')

################
### CONTENTS ###
################

### SETUP ###

	### libraries, functions, and definitions ###
	### define study region: humid forest ###
	### crop forest maps to study region ###

### MODELING DEFORESTATION ###

	### collate data for modeling deforestation ###
	### estimate deforestation VOLUME (vs rate) ###
	### Vielledent model of deforestation RATE ###
	### Vielledent model of deforestation LOCATION ###
	### evaluate deforestation LOCATION models ###
	### project deforestation into future ###
	### create maps of future forest cover assuming no additional loss in protected areas after 2014 ###
	### calculate future forest fragmentation class ###
	### compile forest cover and fragmentation statistics ###
	### create display maps of forest cover ###
	
	### compile GIS files for faster plotting ###
	### create display maps of just Madagascar for presentations ###
	### create display maps of PAST forest cover for presentations ###
	### create display maps of FUTURE forest cover for presentations ###
	### create display maps of forest cover in MOBOT reserves for presentations ###
	
	### create 3D display maps of forest cover ###
	### create display maps of deforestation probability ###
	### create display maps of forest fragmentation class ###

#############################################
### libraries, functions, and definitions ###
#############################################

	memory.limit(memory.limit() * 2^30)
	rm(list=ls())
	gc()

	library(omnibus)
	library(enmSdm)
	library(legendary)
	library(fasterRaster)

	library(sp)
	library(graticule)
	library(raster)
	library(parallel)
	library(dismo)
	library(rgeos)

	# library(phcfM)

	library(purrr)
	library(doBy)
	# library(rayshader)
	library(rgl)
	library(scales)
	library(tictoc)

	say(date())

	###############
	### options ###
	###############

		options(stringsAsFactors=FALSE)
		rasterOptions(format='GTiff', overwrite=TRUE, tmpdir='C:/ecology/!Scratch/_raster')

	#################
	### variables ###
	#################

		madEaProj <- '+init=epsg:32738'

		bioclims <- c(4, 10, 15, 16) # BIOCLIM variables to use
		gcms <- c('CCSM4', 'MIROC-ESM-CHEM', 'GISS-E2-R') # GCMs to use

		longLat <- c('LONG_WGS84', 'LAT_WGS84') # BIOCLIM variables to use

		neighs <- c(3, 5) # neighborhood size (in # of cells) to use to calculate forest loss

		grassDir <- c('C:/OSGeo4W64/', 'grass-7.4.1', 'osgeo4W')
		
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

		# equation for annualized forest loss rate from Puyravaud, J-P. 2003 Standardizing the calculation of the annual rate of deforestation. Forest Ecology and Management 177:593-596.
		# t1, t1	Time 1 and 2
		# A1, A2	Area of forest in times 1 and 2
		compoundInterestLaw <- function(t1, t2, A1, A2) (1 / (t2 - t1)) * log(A2 / A1)

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

# say('#########################################')
# say('### define study region: humid forest ###')
# say('#########################################')

	# say('Ecoregions of Madagascar from citations in Vieilledent.')

	# say('### UTM 38S ###')

		# # get humid forest ecoregion
		# ecoregions_utm38s <- shapefile('./Data/Ecoregions of Madagascar - Citations in Vieilledent et al 2018 Fig 1/madagascar_ecoregion_tenaizy_38s')
		# ecoregions_utm38s <- gBuffer(ecoregions_utm38s, width=0, byid=TRUE) # fixed problem with self-intersection in shapefile

		# humidForest_utm38s <- ecoregions_utm38s[ecoregions_utm38s$nature == 'Ecoregion Est', ]
		# humidForest_utm38s <- sp::spTransform(humidForest_utm38s, CRS(madEaProj))

		# # remove holes
		# p4s <- projection(humidForest_utm38s)
		# humidNoHole <- SpatialPolygons(list(Polygons(list(humidForest_utm38s@polygons[[1]]@Polygons[[1]]), ID=1)))
		# humidForest_utm38s <- SpatialPolygonsDataFrame(humidNoHole, data=humidForest_utm38s@data, match.ID=FALSE)
		# projection(humidForest_utm38s) <- p4s

		# # remove disjunct areas in north
		# humidForest_utm38s <- multiPart(humidForest_utm38s)
		# humidForest_utm38s$area_km2 <- gArea(humidForest_utm38s, byid=TRUE) / (1000^2)
		# humidForest_utm38s <- humidForest_utm38s[which.max(humidForest_utm38s$area_km2), ]

		# # buffer to add all presences/absences for Varecia
		# varecia <- read.csv('./Data/Varecia/02 Master__Varecia_Point_Data_2017 - Removed and Corrected Badly Georeferenced Sites and Zahamena site by EEL on 3-10-2000.csv')
		# varecia <- SpatialPointsDataFrame(varecia[ , longLat], data=varecia, proj4string=getCRS('wgs84', TRUE))

		# # get distance between outlying site and ecoregion
		# varecia_utm38s <- sp::spTransform(varecia, CRS(madEaProj))
		# distToHumid <- gDistance(varecia_utm38s, humidForest_utm38s, byid=TRUE)

		# maxDist <- max(distToHumid)

		# # add buffer to humid forest ecoregion and crop
		# bufferSize <- 3 * maxDist
		# humidForestBuffer_utm38s <- gBuffer(humidForest_utm38s, width=bufferSize)

		# humidForestBuffer_utm38s <- crop(humidForestBuffer_utm38s, ecoregions_utm38s)

		# # create mask
		# forest <- raster('./Data/Forest - Vieilledent et al 2018/ORIGINALS/for2000.tif')
		# proj4string(forest) <- madEaProj
		# forest <- crop(forest, humidForestBuffer_utm38s)
		# humidForestBufferMask_utm38s <- fasterRasterize(humidForestBuffer_utm38s, forest, grassDir=grassDir)
		# names(humidForestBufferMask_utm38s) <- 'humidForestBufferMask_utm38s'

		# # country polygon
		# madagascar_wgs84 <- getData('GADM', country='MDG', level=0, path='./Study Region & Masks')
		# madagascar_wgs84 <- gadm[gadm$NAME_0 == 'Madagascar', ]
		# madagascar_utm38s <- sp::spTransform(madagascar_wgs84, CRS(madEaProj))


		# # Jesus saves
		# dirCreate('./Study Region & Masks/UTM 38S 30-m Resolution')

		# save(humidForest_utm38s, file='./Study Region & Masks/UTM 38S 30-m Resolution/Eastern Humid Forest Polygon.RData')
		# save(humidForestBuffer_utm38s, file='./Study Region & Masks/UTM 38S 30-m Resolution/Eastern Humid Forest Polygon Buffer.RData')

		# writeRaster(humidForestBufferMask_utm38s, './Study Region & Masks/UTM 38S 30-m Resolution/humidForestBufferMask_utm38s', datatype='INT1U')

		# save(madagascar_utm38s, file='./Study Region & Masks/UTM 38S 30-m Resolution/Madagascar from GADM.RData')

	# say('### WGS 84 30 arcsec resolution ###')

		# # project to WGS84
		# humidForest_wgs84 <- sp::spTransform(humidForest_utm38s, getCRS('wgs84', TRUE))
		# humidForestBuffer_wgs84 <- sp::spTransform(humidForestBuffer_utm38s, getCRS('wgs84', TRUE))

		# # create mask
		# wc <- raster('D:/Ecology/Climate/WORLDCLIM Ver 2 Rel June 1 2016/30 arcsec 1970 to 2000/wc2.0_bio_30s_01.tif')
		# wc <- crop(wc, humidForestBuffer_wgs84)
		# humidForestBufferMask_wgs84 <- rasterize(humidForestBuffer_wgs84, wc)
		# humidForestBufferMask_wgs84 <- 1 + 0 * humidForestBufferMask_wgs84
		# names(humidForestBufferMask_wgs84) <- 'humidForestBufferMask_wgs84'

		# # Jesus saves
		# dirCreate('./Study Region & Masks/WGS84 30-arcsec Resolution')

		# save(humidForest_wgs84, file='./Study Region & Masks/WGS84 30-arcsec Resolution/Eastern Humid Forest Polygon.RData')
		# save(humidForestBuffer_wgs84, file='./Study Region & Masks/WGS84 30-arcsec Resolution/Eastern Humid Forest Polygon Buffer.RData')

		# writeRaster(humidForestBufferMask_wgs84, './Study Region & Masks/WGS84 30-arcsec Resolution/humidForestBufferMask_wgs84', datatype='INT1U')

		# save(madagascar_wgs84, file='./Study Region & Masks/WGS84 30-arcsec Resolution/Madagascar from GADM.RData')

# say('########################################')
# say('### crop forest maps to study region ###')
# say('########################################')

	# years <- c(2000, 2005, 2010, 2014)

	# humidForestBufferMask_utm38s <- raster('./Study Region & Masks/UTM 38S 30-m Resolution/humidForestBufferMask_utm38s.tif')

	# for (year in years) {

		# say(year)
		# rast <- raster(paste0('./Data/Forest - Vieilledent et al 2018/ORIGINALS/for', year, '.tif'))
		# rast <- crop(rast, humidForestBufferMask_utm38s)
		# rast <- rast * humidForestBufferMask_utm38s
		# names(rast) <- paste0('forest', year)

		# writeRaster(rast, paste0('./Data/Forest - Vieilledent et al 2018/forest', year), datatype='INT1U')
		# rm(rast); gc()

	# }

# say('###############################################')
# say('### collate data for modeling deforestation ###')
# say('###############################################')

	# humidForestBufferMask_utm38s <- raster('./Study Region & Masks/UTM 38S 30-m Resolution/humidForestBufferMask_utm38s.tif')

	# say('### long/lat rasters')

		# ll <- longLatRasters(humidForestBufferMask_utm38s)
		# ll <- ll * humidForestBufferMask_utm38s
		# names(ll) <- c('longitude', 'latitude')

		# dirCreate('./Data/Longitude & Latitude Rasters')
		# writeRaster(ll, './Data/Longitude & Latitude Rasters/longLat_utm38s', datatype='INT4S')

	# say('### distance to nearest settlement')

		# places <- read.csv('./Data/Geographic Places (NIMA from DIVA-GIS)/00 Places.csv')
		# populatedPlaces <- places[places$F_DESIG == 'PPL', ]
		# populatedPlaces <- SpatialPointsDataFrame(cbind(populatedPlaces$LONG, populatedPlaces$LAT), data=populatedPlaces, proj4string=getCRS('wgs84', TRUE))

		# populatedPlaces <- sp::spTransform(populatedPlaces, CRS(madEaProj))

		# rast <- fasterVectToRastDistance(rast=humidForestBufferMask_utm38s, vect=populatedPlaces, grassDir=grassDir)

		# rast <- round(rast)
		# rast <- humidForestBufferMask_utm38s * rast
		# names(rast) <- 'distToSettlementMeters_utm38s'
		# writeRaster(rast, './Data/Geographic Places (NIMA from DIVA-GIS)/distToSettlementMeters_utm38s', datatype='INT4S')

	# say('### distance to major roads')

		# vect <- shapefile('./Data/Roads/MDG_roads')
		# vect <- sp::spTransform(vect, CRS(madEaProj))

		# out <- fasterVectToRastDistance(vect=vect, rast=humidForestBufferMask_utm38s, metric='euclidean', meters=TRUE, grassDir=grassDir)

		# out <- out * humidForestBufferMask_utm38s
		# out <- round(out)
		# names(out) <- 'distToNearestRoadMeters_utm38s'

		# writeRaster(out, './Data/Roads/distToNearestRoadMeters_utm38s', datatype='INT4S')

	# say('### distance to any water body (marine/inland)')

		# ### lakes

			# say('lakes')

			# vect <- shapefile('./Data/Water Bodies/MDG_water_areas_dcw')
			# vect <- sp::spTransform(vect, CRS(projection(madEaProj)))
			# vect <- crop(vect, humidForestBufferMask_utm38s)

			# out <- fasterVectToRastDistance(rast=humidForestBufferMask_utm38s, vect=vect, metric='euclidean', meters=TRUE, grassDir=grassDir)

			# out <- out * humidForestBufferMask_utm38s
			# out <- round(out)
			# names(out) <- 'distToNearestLakeMeters_utm38s'

			# writeRaster(out, './Data/Water Bodies/distToNearestLakeMeters_utm38s', datatype='INT4S')

		# ### rivers

			# say('rivers')

			# vect <- shapefile('./Data/Water Bodies/MDG_water_lines_dcw')
			# vect <- sp::spTransform(vect, CRS(projection(madEaProj)))
			# vect <- crop(vect, humidForestBufferMask_utm38s)

			# out <- fasterVectToRastDistance(rast=humidForestBufferMask_utm38s, vect=vect, metric='euclidean', meters=TRUE, grassDir=grassDir)

			# out <- out * humidForestBufferMask_utm38s
			# out <- round(out)
			# names(out) <- 'distToNearestRiverMeters_utm38s'

			# writeRaster(out, './Data/Water Bodies/distToNearestRiverMeters_utm38s', datatype='INT4S')

		# say('### distance to any inland water ###')

			# water <- raster::stack(c(
				# './Data/Water Bodies/distToNearestLakeMeters_utm38s.tif',
				# './Data/Water Bodies/distToNearestRiverMeters_utm38s.tif'
			# ))

			# beginCluster(7)
				# out <- calc(water, min)
			# endCluster()

			# names(out) <- 'distToNearestInlandWaterMeters_utm38s'
			# writeRaster(out, './Data/Water Bodies/distToNearestInlandWaterMeters_utm38s', datatype='INT4S')

		# say('### distance to coast')

			# load('./Study Region & Masks/UTM 38S 30-m Resolution/Madagascar from GADM.RData')
			# madMask_utm38s <- fasterRasterize(vect=madagascar_utm38s, rast=humidForestBufferMask_utm38s, grassDir=grassDir)

			# # calculate distance to nearest coast
			# out <- fasterRastDistance(madMask_utm38s, metric='euclidean', meters=TRUE, fillNAs=FALSE, grassDir=grassDir)

			# out <- out * humidForestBufferMask_utm38s
			# out <- round(out)
			# names(out) <- 'distToNearestCoastMeters_utm38s'

			# writeRaster(out, './Data/Water Bodies/distToNearestCoastMeters_utm38s', datatype='INT4S')

		# say('### distance to any water')

			# distToWaterBody <- raster::stack(c(
				# './Data/Water Bodies/distToNearestLakeMeters_utm38s.tif',
				# './Data/Water Bodies/distToNearestRiverMeters_utm38s.tif',
				# './Data/Water Bodies/distToNearestCoastMeters_utm38s.tif')
			# )

			# out <- min(distToWaterBody)
			# names(out) <- 'distToNearestWaterMeters_utm38s'

			# writeRaster(out, './Data/Water Bodies/distToNearestWaterBodyMeters_utm38s')

	# say('### topography from GMTED2010')

		# ### elevation

			# elev <- fasterProjectRaster('D:/Ecology/Topography/GMTED2010/7pt5_arcsec/30S030E_20101117_gmted_mea075.tif', template=humidForestBufferMask_utm38s, method='bilinear', grassDir=grassDir)

		# ### slope

			# slope <- fasterTerrain(elev, slope=TRUE, slopeUnits='degrees', grassDir=grassDir)
			# slope <- slope * humidForestBufferMask_utm38s
			# slope <- round(slope)
			# names(slope) <- 'slope'

			# writeRaster(slope, './Data/Topography - GMTED2010/slopeGmted2010_utm38s', dataType='INT1U')

		# ### elevation (crop & save)

			# elev <- elev * humidForestBufferMask_utm38s
			# elev <- round(elev)
			# names(elev) <- 'elevation_gmted2010'

			# dirCreate('./Data/Topography - GMTED2010')
			# writeRaster(elev, './Data/Topography - GMTED2010/elevationGmted2010_utm38s', dataType='INT2S')

	# say('### forest layers ###')

		# forest2000 <- raster('./Data/Forest - Vieilledent et al 2018/forest2000.tif')
		# forest2005 <- raster('./Data/Forest - Vieilledent et al 2018/forest2005.tif')
		# forest2010 <- raster('./Data/Forest - Vieilledent et al 2018/forest2010.tif')
		# forest2014 <- raster('./Data/Forest - Vieilledent et al 2018/forest2014.tif')

		# say('deforestation')

			# # convert NAs to 0
			# beginCluster(7)

				# zeroFunct <- function(x) ifelse(is.na(x), 0, x)

				# forest2000zeros <- clusterR(forest2000, calc, args=list(fun=zeroFunct))
				# forest2005zeros <- clusterR(forest2005, calc, args=list(fun=zeroFunct))
				# forest2010zeros <- clusterR(forest2010, calc, args=list(fun=zeroFunct))
				# forest2014zeros <- clusterR(forest2014, calc, args=list(fun=zeroFunct))

			# endCluster()

			# forest2000zeros <- humidForestBufferMask_utm38s * forest2000zeros
			# forest2005zeros <- humidForestBufferMask_utm38s * forest2005zeros
			# forest2010zeros <- humidForestBufferMask_utm38s * forest2010zeros
			# forest2014zeros <- humidForestBufferMask_utm38s * forest2014zeros

			# deforest2005 <- forest2000 - forest2005zeros
			# deforest2010 <- forest2005 - forest2010zeros
			# deforest2014 <- forest2010 - forest2014zeros

			# # convert 0s (deforested) to NAs
			# beginCluster(7)

				# fx <- function(x) ifelse(x == 1, 1, NA)
				# deforest2005 <- clusterR(deforest2005, calc, args=list(fun=fx))
				# deforest2010 <- clusterR(deforest2010, calc, args=list(fun=fx))
				# deforest2014 <- clusterR(deforest2014, calc, args=list(fun=fx))

			# endCluster()

			# names(deforest2005) <- 'deforest2005_utm38s'
			# names(deforest2010) <- 'deforest2010_utm38s'
			# names(deforest2014) <- 'deforest2014_utm38s'

			# writeRaster(deforest2005, './Data/Forest - Vieilledent et al 2018/deforest2005_utm38s', datatype='INT1U')
			# writeRaster(deforest2010, './Data/Forest - Vieilledent et al 2018/deforest2010_utm38s', datatype='INT1U')
			# writeRaster(deforest2014, './Data/Forest - Vieilledent et al 2018/deforest2014_utm38s', datatype='INT1U')

		# say('### distance to deforested patch ###')

			# say('2005')
			# distRast <- fasterRastDistance(deforest2005, metric='euclidean', meters=TRUE, fillNAs=TRUE, grassDir=grassDir)
			# distRast <- distRast * humidForestBufferMask_utm38s
			# distRast <- round(distRast)
			# names(distRast) <- 'distToDefoMeters2005_utm38s'
			# writeRaster(distRast, './Data/Forest - Vieilledent et al 2018/distToDefoMeters2005_utm38s', datatype='INT4S')

			# say('2010')
			# distRast <- fasterRastDistance(deforest2010, metric='euclidean', meters=TRUE, fillNAs=TRUE, grassDir=grassDir)
			# distRast <- distRast * humidForestBufferMask_utm38s
			# distRast <- round(distRast)
			# names(distRast) <- 'distToDefoMeters2010_utm38s'
			# writeRaster(distRast, './Data/Forest - Vieilledent et al 2018/distToDefoMeters2010_utm38s', datatype='INT4S')

			# say('2014')
			# distRast <- fasterRastDistance(deforest2014, metric='euclidean', meters=TRUE, fillNAs=TRUE, grassDir=grassDir)
			# distRast <- distRast * humidForestBufferMask_utm38s
			# distRast <- round(distRast)
			# names(distRast) <- 'distToDefoMeters2014_utm38s'
			# writeRaster(distRast, './Data/Forest - Vieilledent et al 2018/distToDefoMeters2014_utm38s', datatype='INT4S')

		# say('### distance to forest edge ###')

			# say('2000')
			# distRast <- fasterRastDistance(forest2000, metric='euclidean', meters=TRUE, fillNAs = FALSE, grassDir=grassDir)
			# distRast <- round(distRast)
			# distRast <- distRast * humidForestBufferMask_utm38s
			# names(distRast) <- 'distToForestEdgeMeters2000_utm38s'
			# writeRaster(distRast, './Data/Forest - Vieilledent et al 2018/distToForestEdgeMeters2000_utm38s', datatype='INT4S')

			# say('2005')
			# distRast <- fasterRastDistance(forest2005, metric='euclidean', meters=TRUE, fillNAs=FALSE, grassDir=grassDir)
			# distRast <- round(distRast)
			# distRast <- distRast * humidForestBufferMask_utm38s
			# names(distRast) <- 'distToForestEdgeMeters2005_utm38s'
			# writeRaster(distRast, './Data/Forest - Vieilledent et al 2018/distToForestEdgeMeters2005_utm38s', datatype='INT4S')

			# say('2010')
			# distRast <- fasterRastDistance(forest2010, metric='euclidean', meters=TRUE, fillNAs=FALSE, grassDir=grassDir)
			# distRast <- round(distRast)
			# distRast <- distRast * humidForestBufferMask_utm38s
			# names(distRast) <- 'distToForestEdgeMeters2010_utm38s'
			# writeRaster(distRast, './Data/Forest - Vieilledent et al 2018/distToForestEdgeMeters2010_utm38s', datatype='INT4S')

			# say('2014')
			# distRast <- fasterRastDistance(forest2014, metric='euclidean', meters=TRUE, fillNAs=FALSE, grassDir=grassDir)
			# distRast <- round(distRast)
			# distRast <- distRast * humidForestBufferMask_utm38s
			# names(distRast) <- 'distToForestEdgeMeters2014_utm38s'
			# writeRaster(distRast, './Data/Forest - Vieilledent et al 2018/distToForestEdgeMeters2014_utm38s', datatype='INT4S')

		# say('### forest fragmentation ###')

			# cores <- 4

			# say('2000')

				# frag2000 <- fasterFragmentation(rast = forest2000zeros, size = 5, pad = TRUE, padValue = NA, calcDensity = TRUE, calcConnect = TRUE, calcClass = TRUE, na.rm = TRUE, undet = 'perforated', cores = cores, forceMulti = TRUE)

				# frag2000 <- frag2000 * humidForestBufferMask_utm38s
				# frag2000[['density']] <- round(100 * frag2000[['density']])
				# frag2000[['connect']] <- round(100 * frag2000[['connect']])

				# names(frag2000) <- c('forestFragClass2000_utm38s', 'forestFragDensity2000_utm38s', 'forestFragConnect2000_utm38s')

				# writeRaster(frag2000[['forestFragClass2000_utm38s']], './Data/Forest - Vieilledent et al 2018/forestFragClass2000_utm38s', datatype='INT1U')
				# writeRaster(frag2000[['forestFragDensity2000_utm38s']], './Data/Forest - Vieilledent et al 2018/forestFragDensity2000_utm38s', datatype='INT1U')
				# writeRaster(frag2000[['forestFragConnect2000_utm38s']], './Data/Forest - Vieilledent et al 2018/forestFragConnect2000_utm38s', datatype='INT1U')

			# say('2005')

				# frag2005 <- fasterFragmentation(rast = forest2005zeros, size = 5, pad = TRUE, padValue = NA, calcDensity = TRUE, calcConnect = TRUE, calcClass = TRUE, na.rm = TRUE, undet = 'perforated', cores = cores, forceMulti = TRUE)

				# frag2005[['density']] <- round(100 * frag2005[['density']])
				# frag2005[['connect']] <- round(100 * frag2005[['connect']])
				# frag2005 <- frag2005 * humidForestBufferMask_utm38s

				# names(frag2005) <- c('forestFragClass2005_utm38s', 'forestFragDensity2005_utm38s', 'forestFragConnect2005_utm38s')

				# writeRaster(frag2005[['forestFragClass2005_utm38s']], './Data/Forest - Vieilledent et al 2018/forestFragClass2005_utm38s', datatype='INT1U')
				# writeRaster(frag2005[['forestFragDensity2005_utm38s']], './Data/Forest - Vieilledent et al 2018/forestFragDensity2005_utm38s', datatype='INT1U')
				# writeRaster(frag2005[['forestFragConnect2005_utm38s']], './Data/Forest - Vieilledent et al 2018/forestFragConnect2005_utm38s', datatype='INT1U')

			# say('2010')

				# frag2010 <- fasterFragmentation(rast = forest2010zeros, size = 5, pad = TRUE, padValue = NA, calcDensity = TRUE, calcConnect = TRUE, calcClass = TRUE, na.rm = TRUE, undet = 'perforated', cores = cores, forceMulti = TRUE)

				# frag2010[['density']] <- round(100 * frag2010[['density']])
				# frag2010[['connect']] <- round(100 * frag2010[['connect']])
				# frag2010 <- frag2010 * humidForestBufferMask_utm38s

				# names(frag2010) <- c('forestFragClass2010_utm38s', 'forestFragDensity2010_utm38s', 'forestFragConnect2010_utm38s')

				# writeRaster(frag2010[['forestFragClass2010_utm38s']], './Data/Forest - Vieilledent et al 2018/forestFragClass2010_utm38s', datatype='INT1U')
				# writeRaster(frag2010[['forestFragDensity2010_utm38s']], './Data/Forest - Vieilledent et al 2018/forestFragDensity2010_utm38s', datatype='INT1U')
				# writeRaster(frag2010[['forestFragConnect2010_utm38s']], './Data/Forest - Vieilledent et al 2018/forestFragConnect2010_utm38s', datatype='INT1U')

			# say('2014')

				# frag2014 <- fasterFragmentation(rast = forest2014zeros, size = 5, pad = TRUE, padValue = NA, calcDensity = TRUE, calcConnect = TRUE, calcClass = TRUE, na.rm = TRUE, undet = 'perforated', cores = cores, forceMulti = TRUE)

				# frag2014[['density']] <- round(100 * frag2014[['density']])
				# frag2014[['connect']] <- round(100 * frag2014[['connect']])
				# frag2014 <- frag2014 * humidForestBufferMask_utm38s

				# names(frag2014) <- c('forestFragClass2014_utm38s', 'forestFragDensity2014_utm38s', 'forestFragConnect2014_utm38s')

				# writeRaster(frag2014[['forestFragClass2014_utm38s']], './Data/Forest - Vieilledent et al 2018/forestFragClass2014_utm38s', datatype='INT1U')
				# writeRaster(frag2014[['forestFragDensity2014_utm38s']], './Data/Forest - Vieilledent et al 2018/forestFragDensity2014_utm38s', datatype='INT1U')
				# writeRaster(frag2014[['forestFragConnect2014_utm38s']], './Data/Forest - Vieilledent et al 2018/forestFragConnect2014_utm38s', datatype='INT1U')

	# say('### protected areas ###')

		# pa <- shapefile('./Data/Protected Areas/WDPA_Sept2018_MDG-shapefile-polygons')

		# # rasterize
		# pa_utm38s <- sp::spTransform(pa, CRS(madEaProj))

		# # crop to humid forest
		# humidForestBufferMask_utm38s <- raster('./Study Region & Masks/UTM 38S 30-m Resolution/humidForestBufferMask_utm38s.tif')

		# pa_utm38s <- crop(pa_utm38s, humidForestBufferMask_utm38s)

		# paRast_utm38s <- fasterRasterize(pa_utm38s, humidForestBufferMask_utm38s, grassDir=grassDir)

		# paRast_utm38s <- humidForestBufferMask_utm38s * paRast_utm38s
		# names(paRast_utm38s) <- 'protectedAreas_utm38s'

		# writeRaster(paRast_utm38s, './Data/Protected Areas/wdpa_utm38s')

# say('##################################################')
# say('### define areas of deforestation for modeling ###')
# say('##################################################')

	# say('Divide forest/deforestation into geographic clusters using a clustering algorithm. Model each separately.')

	# humidForestBufferMask_utm38s <- raster('./Study Region & Masks/UTM 38S 30-m Resolution/humidForestBufferMask_utm38s.tif')

	# # generate sample of points from areas deforested 2000-2014

		# forest2000 <- raster('./Data/Forest - Vieilledent et al 2018/forest2000.tif')

		# forestSamples <- sampleRast(forest2000, n=10000, adjArea=FALSE, replace=FALSE, prob=FALSE)

	# # cluster deforested sites

		# forestSamples <- SpatialPoints(forestSamples, CRS(madEaProj))

		# forestSamples <- sp::spTransform(forestSamples, getCRS('wgs84', TRUE))
		# dists <- pointDist(forestSamples)
		# dists <- as.dist(dists)

		# library(cluster)
		# groups <- agnes(dists)

		# png('./Deforestation Models/Dendrogram of Forest Clusters Calculated Using UPGMA.png', width=900, height=900)
			# plot(groups, main='Forest Clusters')
		# dev.off()

		# numFolds <- 6
		# gFolds <- geoFold(forestSamples, numFolds)

	# # assign cells to each group

		# dirCreate('./Deforestation Models')

		# for (i in 1:numFolds) {

			# say('Calculating distance to fold ', i, '...')

			# foldSamples <- defoSamples[gFolds == i]
			# center <- gCentroid(foldSamples)

			# distToFold <- fasterVectToRastDistance(humidForestBufferMask_utm38s, center)
			# distToFold <- round(distToFold)
			# names(distToFold) <- paste0('distToFold_', i)
			# writeRaster(distToFold, paste0('./Deforestation Models/distToFold_', i), datatype='INT4S')

		# }

		# distToFolds <- raster::stack(listFiles('./Deforestation Models', pattern='distToFold'))

		# beginCluster(7)
			# foldAssign <- calc(distToFolds, which.min)
		# endCluster()

		# names(foldAssign) <- 'deforestationRegion'
		# writeRaster(foldAssign, './Deforestation Models/deforestationRegion', datatype='INT1U')


# say('############################################################')
# say('### generate samples for modeling deforestation LOCATION ###')
# say('############################################################')

	# say('Predictors:')
	# say('   elevation')
	# say('   slope')
	# say('   forest fragmentation index')
	# say('   distance to forest edge')
	# say('   distance to previously deforested patch')
	# say('   distance to road')
	# say('   distance to settlement')
	# say('   distance to nearest inland water body')
	# say('   distance to nearest coast')
	# say('   un/protected')
	# say('   longitude and latitude')
	# say('   human population density')

	# # number of forested sites in starting time period... getting a little extra to account for NAs
	# numTrainSites <- 40000
	# numTestSites <- 10000
	# numSites <- numTrainSites + numTestSites

	# periodStarts <- c(2000, 2005, 2010)
	# periodEnds <- c(2005, 2010, 2014)

	# dirCreate('./Deforestation Models')

	# ### by TIME PERIOD
	# for (period in 1:3) {

		# periodStart <- periodStarts[period]
		# periodEnd <- periodEnds[period]

		# say('Extracting deforestation predictors for ', periodStart, ' to ', periodEnd)

		# ### generate sample sites and response variable (hard to place sites on the sparse forest rasters so getting more than needed)

			# forestStart <- raster(paste0('./Data/Forest - Vieilledent et al 2018/forest', periodStart, '.tif'))
			# forestEnd <- raster(paste0('./Data/Forest - Vieilledent et al 2018/forest', periodEnd, '.tif'))

			# sites <- sampleRast(forestStart, 1.05 * numSites, adjArea=FALSE, replace=FALSE, prob=FALSE)

			# # create response variable (0 = not deforested, 1 = deforested)
			# coverStart <- extract(forestStart, sites)
			# coverEnd <- extract(forestEnd, sites)

			# defo <- rep(0, numTrainSites + numTestSites)
			# defo[is.na(coverEnd)] <- 1

		# ### extract predictors

		# preds <- as.data.frame(sites)
		# names(preds) <- c('longitude', 'latitude')

			# # elevation

				# rast <- raster('./Data/Topography - GMTED2010/elevationGmted2010_utm38s.tif')
				# preds <- cbind(preds, extract(rast, sites))
				# names(preds)[ncol(preds)] <- 'elevation'

			# # slope

				# rast <- raster('./Data/Topography - GMTED2010/slopeGmted2010_utm38s.tif')
				# preds <- cbind(preds, extract(rast, sites))
				# names(preds)[ncol(preds)] <- 'slope'

			# # distance to nearest settlement

				# rast <- raster('./Data/Geographic Places (NIMA from DIVA-GIS)/distToSettlementMeters_utm38s.tif')
				# preds <- cbind(preds, extract(rast, sites))
				# names(preds)[ncol(preds)] <- 'distToSettlement'

			# # distance to nearest road

				# rast <- raster('./Data/Roads/distToNearestRoadMeters_utm38s.tif')
				# preds <- cbind(preds, extract(rast, sites))
				# names(preds)[ncol(preds)] <- 'distToRoad'

			# # distance to nearest inland water body

				# rast <- raster('./Data/Water Bodies/distToNearestInlandWaterMeters_utm38s.tif')
				# preds <- cbind(preds, extract(rast, sites))
				# names(preds)[ncol(preds)] <- 'distToInlandWater'

			# # distance to nearest coast

				# rast <- raster('./Data/Water Bodies/distToNearestCoastMeters_utm38s.tif')
				# preds <- cbind(preds, extract(rast, sites))
				# names(preds)[ncol(preds)] <- 'distToCoast'

			# # un/protected

				# rast <- raster('./Data/Protected Areas/wdpa_utm38s.tif')
				# preds <- cbind(preds, extract(rast, sites))
				# names(preds)[ncol(preds)] <- 'protectedArea'
				# preds$protectedArea[is.na(preds$protectedArea)] <- 0

			# # forest fragmentation class

				# rast <- raster(paste0('./Data/Forest - Vieilledent et al 2018/forestFragClass', periodStart, '_utm38s.tif'))
				# preds <- cbind(preds, extract(rast, sites))
				# names(preds)[ncol(preds)] <- 'fragClass'

			# # distance to nearest deforested patch (not included for 2000-2005)

				# if (periodStart != 2000) {

					# rast <- raster(paste0('./Data/Forest - Vieilledent et al 2018/distToDefoMeters', periodStart, '_utm38s.tif'))
					# preds <- cbind(preds, extract(rast, sites))
					# names(preds)[ncol(preds)] <- 'distToDefo'

				# } else {

					# preds$distToDefo <- NA

				# }

			# # distance to nearest forest edge

				# rast <- raster(paste0('./Data/Forest - Vieilledent et al 2018/distToForestEdgeMeters', periodStart, '_utm38s.tif'))
				# preds <- cbind(preds, extract(rast, sites))
				# names(preds)[ncol(preds)] <- 'distToForestEdge'

			# # human population

				# sites_utm38s <- SpatialPoints(sites, CRS(madEaProj))
				# sites_wgs84 <- sp::spTransform(sites_utm38s, getCRS('wgs84', TRUE))

				# rast <- raster(paste0('./Data/Population/popDensity', periodStart, '_wgs84.tif'))
				# preds <- cbind(preds, extract(rast, sites_wgs84))
				# names(preds)[ncol(preds)] <- 'popDensity'

		# # collate data

		# intervals <- data.frame(interval=rep(periodEnd - periodStart, nrow(preds)))

		# thisData <- cbind(defo, intervals, preds)

		# thisData$slope[is.na(thisData$slope)] <- 0

		# nas <- naRows(thisData)
		# if (length(nas) > 0) thisData <- thisData[1:numSites, ]

		# thisData$startEndYear <- paste(periodStart, periodEnd)
		# thisData$period <- period

		# trainTest <- data.frame(trainTest = c(rep('train', numTrainSites), rep('test', numTestSites)))
		# thisData <- as.data.frame(thisData)
		# thisData <- cbind(trainTest, thisData)

		# collatedData <- if (exists('collatedData')) {
			# rbind(collatedData, thisData)
		# } else {
			# thisData
		# }

	# } # next time period

	# collatedData$fragClass[collatedData$fragClass == 5] <- 4

	# trainData <- collatedData[collatedData$trainTest == 'train', ]
	# testData <- collatedData[collatedData$trainTest == 'test', ]

	# save(trainData, file='./Deforestation Models/Training Sites & Predictors Located in Humid Forest UTM 38S.RData')
	# save(testData, file='./Deforestation Models/Test Sites & Predictors Located in Humid Forest UTM 38S.RData')

# say('###############################################')
# say('### estimate deforestation VOLUME (vs rate) ###')
# say('###############################################')

	# say('We are estimating deforestation VOLUME (vs rate) as the total area deforested each year as a function of total population size. Using rate assumes that the same *proportion* of forest is lost each year, whereas a static population should require the same amount of fiber product per year, regardless of stock. We will assume that the total amount of deforestation occurring in a year is a linear function of total population size such that: D = alpha * N0 + alpha * N1 + alpha * N2 + alpha * N3 + alpha * N4 = alpha * SUM N, where D is total area deforested between year Y0 and Y4, alpha is a constant, and N0, N1, N3, ... is the total population in years 0, 1, 2, ... for a given period.', breaks=80)

	# load('./Data/Population Forecasts (UN)/WPP2017_PopulationBySingleAgeSex.RData')
	# madPop <- worldPop[worldPop$Location == 'Madagascar', ]

	# madPop <- summaryBy(PopTotal ~ Time, data=madPop, FUN=sum)
	# names(madPop) <- c('year', 'population')

	# deforest2005 <- raster('./Data/Forest - Vieilledent et al 2018/deforest2005_utm38s.tif')
	# deforest2010 <- raster('./Data/Forest - Vieilledent et al 2018/deforest2010_utm38s.tif')
	# deforest2014 <- raster('./Data/Forest - Vieilledent et al 2018/deforest2014_utm38s.tif')

	# D1 <- cellStats(deforest2005, 'sum')
	# D2 <- cellStats(deforest2010, 'sum')
	# D3 <- cellStats(deforest2014, 'sum')

	# # re-weight first and last years assuming population data was taken in middle of each year
	# N1 <- madPop$population[madPop$year >= 2000 & madPop$year <= 2004]
	# N1[1] <- 0.5 * N1[1]
	# N1[length(N1)] <- 0.5 * N1[length(N1)]

	# N2 <- madPop$population[madPop$year >= 2005 & madPop$year <= 2010]
	# N2[1] <- 0.5 * N2[1]
	# N2[length(N2)] <- 0.5 * N2[length(N2)]

	# N3 <- madPop$population[madPop$year >= 2010 & madPop$year <= 2014]
	# N3[1] <- 0.5 * N3[1]
	# N3[length(N3)] <- 0.5 * N3[length(N3)]

	# alpha1 <- D1 / sum(N1)
	# alpha2 <- D2 / sum(N2)
	# alpha3 <- D3 / sum(N3)

	# ### predict forest amount loss per year... comes to more than total forest remaining in 2014, even if using *minimum* amount rate
	# # MINIMUM
	# defoAmount <- min(alpha1, alpha2, alpha3) * madPop$population[madPop$year >= 2015 & madPop$year <= 2080]
	# defoAmount <- data.frame(year=2015:2080, amount_30x30m=defoAmount)

	# dirCreate('./Figures & Tables/Forest - Deforestation Model Parameterization')
	# write.csv(defoAmount, './Figures & Tables/Forest - Deforestation Model Parameterization/Deforested Amount per Year Using Minimum Loss Rate.csv', row.names=FALSE)

	# ### predict forest amount loss per year... comes to more than total forest remaining in 2014, even if using *minimum* amount rate
	# # MEAN
	# defoAmount <- mean(alpha1, alpha2, alpha3) * madPop$population[madPop$year >= 2015 & madPop$year <= 2080]
	# defoAmount <- data.frame(year=2015:2080, amount_30x30m=defoAmount)

	# dirCreate('./Figures & Tables/Forest - Deforestation Model Parameterization')
	# write.csv(defoAmount, './Figures & Tables/Forest - Deforestation Model Parameterization/Deforested Amount per Year Using Mean Loss Rate.csv', row.names=FALSE)

	# ### predict forest amount loss per year... comes to more than total forest remaining in 2014, even if using *minimum* amount rate
	# # MAXIMUM
	# defoAmount <- max(alpha1, alpha2, alpha3) * madPop$population[madPop$year >= 2015 & madPop$year <= 2080]
	# defoAmount <- data.frame(year=2015:2080, amount_30x30m=defoAmount)

	# dirCreate('./Figures & Tables/Forest - Deforestation Model Parameterization')
	# write.csv(defoAmount, './Figures & Tables/Forest - Deforestation Model Parameterization/Deforested Amount per Year Using Maximum Loss Rate.csv', row.names=FALSE)

	# say('So... even if we use the minimum amount loss rate across the three periods, the total forest lost by 2080 is > all remaining forest.')

# say('##############################################')
# say('### Vielledent model of deforestation RATE ###')
# say('##############################################')

	# load('./Deforestation Models/Training Sites & Predictors Located in Humid Forest UTM 38S.RData')

	# model <- deforestation(defo ~ 1, interval=trainData$interval, data=trainData, burnin=8000, mcmc=2000, thin=1)

	# mcmcMean <- as.matrix(model$mcmc)
	# thetaMean <- mean(inv.logit(mcmcMean[ , 1])) # posterior mean of the annual deforestation rate

	# dirCreate('./Figures & Tables/Forest - Deforestation Model Parameterization')
	# write.csv(thetaMean, './Figures & Tables/Forest - Deforestation Model Parameterization/Mean Interval-Sensitive Deforestation Rate.csv', row.names=FALSE)

# say('##################################################')
# say('### Vielledent model of deforestation LOCATION ###')
# say('##################################################')

	# load('./Deforestation Models/Training Sites & Predictors Located in Humid Forest UTM 38S.RData')

	# # trainData <- trainData[!is.na(trainData$distToDefo), ]

	# # transform
	# trainData$longitude <- stretchMinMax(trainData$longitude, lower=-1, upper=1)
	# trainData$latitude <- stretchMinMax(trainData$latitude, lower=-1, upper=1)
	# trainData$distToSettlement <- log10(trainData$distToSettlement + 0.1)
	# trainData$distToRoad <- log10(trainData$distToRoad + 0.1)
	# trainData$distToCoast <- log10(trainData$distToCoast + 0.1)
	# trainData$distToInlandWater <- log10(trainData$distToInlandWater + 0.1)
	# trainData$distToDefo <- log10(trainData$distToDefo + 1)
	# trainData$distToForestEdge <- log10(trainData$distToForestEdge + 0.1)

	# variables <- c('longitude', 'latitude', 'elevation', 'slope', 'distToSettlement', 'distToRoad', 'distToInlandWater', 'distToCoast', 'as.factor(protectedArea)', 'as.factor(fragClass)', 'distToDefo', 'distToForestEdge', 'popDensity')

	# trainData <- trainData[trainData$period != 1, ]

	# ### evaluate predictor importance

		# say('evaluate predictor importance')

		# # formulas for all variables except one term
		# forms <- list()

		# for (v in seq_along(variables)) {
			# forms[[v]] <- formula(paste('defo ~ ', paste(variables[-v], collapse='+'), sep=''))
		# }

		# # full model and null model
		# forms[[length(forms) + 1]] <- formula(paste('defo ~ ', paste(variables, collapse='+'), sep=''))
		# forms[[length(forms) + 1]] <- formula('defo ~ 1')

		# # deviance
		# dev <- rep(NA, length(forms))

		# # train successive models
		# for (i in seq_along(forms)) {

			# say(forms[[i]])

			# # train model
			# model <- deforestation(forms[[i]], data=trainData, interval=trainData$interval, burnin=8000, mcmc=2000)
			# mcmc <- as.matrix(model$mcmc)
			# gammaHat <- apply(mcmc, 2 , mean)

			# # predict to train data
			# mfFixed <- model.frame(formula=forms[[i]], data=trainData)
			# X <- model.matrix(attr(mfFixed, 'terms'), data=mfFixed)
			# thetaHat <- inv.logit(X %*% gammaHat)
			# thetaPrim <- 1 - (1 - thetaHat)^(trainData$interval)

			# dev[i] <- -2 * sum(dbinom(trainData$defo, 1 , thetaPrim, 1))

		# }

		# # model comparison
		# modelComp <- as.data.frame(matrix(NA, nrow=length(forms) + 1, ncol=2))
		# names(modelComp) <- c('model', 'deviance')
		# modelComp$model <-  c(paste0('-', variables), 'full', 'null', 'saturated')
		# modelComp$deviance <- c(round(dev, 3), 0)

		# dirCreate('./Figures & Tables/Forest - Deforestation Model Parameterization')
		# write.csv(modelComp, './Figures & Tables/Forest - Deforestation Model Parameterization/Deforestation Location Model Comparison.csv', row.names=FALSE)

		# # variable importance (higher deviance ==> greater importance)
		# varImp <- as.data.frame(matrix(NA, nrow=length(forms) - 1, ncol=3))
		# names(varImp) <- c('model', 'deviance', 'varImp')
		# varImp$model <- c(paste0('-', variables), 'full')
		# varImp$deviance <- round(dev, 3)[1:(length(forms) - 1)]
		# varImp$varImp <- (varImp$deviance - varImp$deviance[varImp$model == 'full'])

		# write.csv(varImp, './Figures & Tables/Forest - Deforestation Model Parameterization/Deforestation Location Variable Importance.csv', row.names=FALSE)

	# ### preliminary final model -- used to explore selected variables

		# say('preliminary final model with selected variables')

		# variables <- varImp$model[varImp$varImp > 0]
		# variables <- substr(variables, 2, nchar(variables))

		# form <- paste0('defo ~ ', paste0(variables, collapse=' + '))

		# model <- deforestation(form, interval=trainData$interval, data=trainData, burnin=8000, mcmc=2000, verbose=1)

		# summary(model$mcmc)
		# plot(model$mcmc)

	# ## full model with selected variables

		# say('full model with selected variables')

		# # removing distToCoast, distToInlandWater bc unexpected relationships to deforestation

		# # removing fragmentation class bc correlated with distToForestEdge and because fragClass has non-zero coefficients only for frag class 4 (positive relationship)

		# # removing popDensity bc has zero coeff
		# variables <- c('longitude', 'elevation', 'slope', 'as.factor(protectedArea)', 'distToDefo', 'distToForestEdge')

		# form <- paste0('defo ~ ', paste0(variables, collapse=' + '))

		# model <- deforestation(form, interval=trainData$interval, data=trainData, burnin=8000, mcmc=2000, verbose=1)
		# save(model, file='./Deforestation Models/Vielledent Variable Interval GLM on Deforestation - Final Model.RData')

		# sink('./Figures & Tables/Forest - Deforestation Model Parameterization/Vielledent Variable Interval GLM on Deforestation - Final Model.txt', split=TRUE)
			# summary(model$mcmc)
			# say('Deviance: ', model$deviance)
		# sink()

		# plot(model$mcmc)

# say('##############################################')
# say('### evaluate deforestation LOCATION models ###')
# say('##############################################')

	# load('./Deforestation Models/Vielledent Variable Interval GLM on Deforestation - Final Model.RData')
	# load('./Deforestation Models/Training Sites & Predictors Located in Humid Forest UTM 38S.RData')
	# load('./Deforestation Models/Test Sites & Predictors Located in Humid Forest UTM 38S.RData')

	# # remove 1st period bc no variable "distToDefo"
	# testData <- testData[!is.na(testData$distToDefo), ]

	# # transform
	# testData$longitude <- 2 * (testData$longitude - min(trainData$longitude)) / (max(trainData$longitude) - min(trainData$longitude)) - 1
	# testData$latitude <- 2 * (testData$latitude - min(trainData$latitude)) / (max(trainData$latitude) - min(trainData$latitude)) - 1
	# testData$distToSettlement <- log10(testData$distToSettlement + 0.1)
	# testData$distToRoad <- log10(testData$distToRoad + 0.1)
	# testData$distToCoast <- log10(testData$distToCoast + 0.1)
	# testData$distToInlandWater <- log10(testData$distToInlandWater + 0.1)
	# testData$distToDefo <- log10(testData$distToDefo + 1)
	# testData$distToForestEdge <- log10(testData$distToForestEdge + 0.1)

	# variables <- c('longitude', 'elevation', 'slope', 'as.factor(protectedArea)', 'distToDefo', 'distToForestEdge')
	# form <- paste0('defo ~ ', paste0(variables, collapse=' + '))

	# betaHats <- apply(model$mcmc, 2, mean)
	# mf <- model.frame(formula=form, data=testData)
	# X <- model.matrix(attr(mf, 'terms'), data=mf)
	# preds <- inv.logit(X %*% betaHats)
	# predsPrime <- c(1 - (1 - preds)^(testData$interval))

	# pres <- preds[trainData$defo == 1]
	# abs <- preds[trainData$defo == 0]

	# evaluate(pres, abs)

	# dirCreate('./Figures & Tables/Forest - Deforestation Model Parameterization')
	# write.csv(auc, './Figures & Tables/Forest - Deforestation Model Parameterization/Deforestation Location Model Performance - AUC.csv', row.names=FALSE)

# say('#########################################')
# say('### project deforestation into future ###')
# say('#########################################')

	# n <- 6
	# say('Starting cluster with ', n, ' cores...')
	# beginCluster(n)

	# ### define basic functions
	# #################################

		# ### log function
		# log10Fx <- function(x, offset=0.1) log10(x + offset)

		# ### convert NAs to 0s
		# naToZeroFx <- function(x) ifelse(is.na(x), 0, x)

		# ### convert 0s to NAs
		# zeroToNaFx <- function(x) ifelse(x == 0, NA, x)

		# ### convert 0s to 1s
		# nonOnesToNaFx <- function(x) ifelse(x == 1, 1, NA)

		# ### predict function
		# # x is predictor raster, beta is coefficient
		# predictFx <- function(x, beta) x * beta

		# ### function to combine making prediction for focal variable and adding it to the prediction raster
		# predictAndSumFx <- function(thisPred, rast, beta) {

			# rast <- clusterR(rast, predictFx, args=list(beta=beta))
			# stacked <- stack(thisPred, rast)
			# sumFx <- function(x) sum(x, na.rm=FALSE)
			# thisPred <- clusterR(stacked, fun=sumFx)
			# thisPred

		# }

	# ### deforestation rates
	# #######################

		# # rates
		# defoRates <- read.csv('./Figures & Tables/Forest - Deforestation Model Parameterization/Mean Interval-Sensitive Deforestation Rate.csv')
		# r <- defoRates$x

		# # amounts
		# amounts <- read.csv('./Figures & Tables/Forest - Deforestation Model Parameterization/Deforested Amount per Year Using Minimum Loss Rate.csv')

	# ### load model
	# ###############

		# load('./Deforestation Models/Vielledent Variable Interval GLM on Deforestation - Final Model.RData')

		# betaHats <- apply(model$mcmc, 2, mean)
		# variables <- names(betaHats)
		# variables <- substr(variables, 6, nchar(variables))
		# names(betaHats) <- variables

	# #####################################################################################
	# ### if predictions using static predictors has already been made, use that raster ###
	# #####################################################################################

	# if (file.exists('./Deforestation Models/Deforestation Location Probability 2015-2080/!defoLocationStaticPredictors_utm38s.tif')) {

		# staticPred <- raster('./Deforestation Models/Deforestation Location Probability 2015-2080/!defoLocationStaticPredictors_utm38s.tif')

	# } else {

		# ### predict to static layers
		# ############################

			# say('Predicting using static predictors:', post=0)

			# humidForestBufferMask_utm38s <- raster('./Study Region & Masks/UTM 38S 30-m Resolution/humidForestBufferMask_utm38s.tif')

			# # intercept
			# say('| intercept', post=0)
			# beta <- betaHats[['(Intercept)']]
			# staticPred <- clusterR(humidForestBufferMask_utm38s, predictFx, args=list(beta=beta))

			# # longitude
			# if ('longitude' %in% variables) {

				# say('longitude', post=0)

				# beta <- betaHats[['longitude']]

				# rast <- raster::stack('./Data/Longitude & Latitude Rasters/longLat_utm38s.tif')
				# rast <- rast[[1]]

				# # rescale
				# load('./Deforestation Models/Training Sites & Predictors Located in Humid Forest UTM 38S.RData')
				# minimum <- min(trainData$longitude)
				# maximum <- max(trainData$longitude)

				# fx <- function(x, minimum, maximum) 2 * (x - minimum) / (maximum - minimum) - 1
				# rast <- clusterR(rast, fun=fx, args=list(minimum=minimum, maximum=maximum))

				# # predict
				# staticPred <- predictAndSumFx(thisPred=staticPred, rast=rast, beta=beta)

			# }

			# # latitude
			# if ('latitude' %in% variables) {

				# say('| latitude', post=0)

				# beta <- betaHats[['latitude']]

				# rast <- raster::stack('./Data/Longitude & Latitude Rasters/longLat_utm38s.tif')
				# rast <- rast[[2]]

				# # rescale
				# load('./Deforestation Models/Training Sites & Predictors Located in Humid Forest UTM 38S.RData')
				# minimum <- min(trainData$latitude)
				# maximum <- max(trainData$latitude)

				# fx <- function(x, minimum, maximum) 2 * (x - minimum) / (maximum - minimum) - 1
				# rast <- clusterR(rast, fun=fx, args=list(minimum=minimum, maximum=maximum))

				# # predict
				# staticPred <- predictAndSumFx(thisPred=staticPred, rast=rast, beta=beta)

			# }

			# # elevation
			# if ('elevation' %in% variables) {

				# say('| elevation', post=0)

				# beta <- betaHats[['elevation']]

				# rast <- raster('./Data/Topography - GMTED2010/elevationGmted2010_utm38s.tif')
				# staticPred <- predictAndSumFx(thisPred=staticPred, rast=rast, beta=beta)

			# }

			# # slope
			# if ('slope' %in% variables) {

				# say('| slope', post=0)

				# beta <- betaHats[['slope']]

				# rast <- raster('./Data/Topography - GMTED2010/slopeGmted2010_utm38s.tif')
				# staticPred <- predictAndSumFx(thisPred=staticPred, rast=rast, beta=beta)

			# }

			# # distance to nearest settlement
			# if ('distToSettlement' %in% variables) {

				# say('| distToSettlement', post=0)

				# beta <- betaHats[['distToSettlement']]

				# rast <- raster('./Data/Geographic Places/distToSettlementMeters_utm38s.tif')
				# rast <- clusterR(rast, log10Fx)
				# staticPred <- predictAndSumFx(thisPred=staticPred, rast=rast, beta=beta)

			# }

			# # distance to nearest road
			# if ('distToRoad' %in% variables) {

				# say('| distToRoad', post=0)

				# beta <- betaHats[['distToRoad']]

				# rast <- raster('./Data/Roads/distToNearestRoadMeters_utm38s.tif')
				# rast <- clusterR(rast, log10Fx)
				# staticPred <- predictAndSumFx(thisPred=staticPred, rast=rast, beta=beta)

			# }

			# # distance to nearest coast
			# if ('distToCoast' %in% variables) {

				# say('| distToCoast', post=0)

				# beta <- betaHats[['distToCoast']]

				# rast <- raster('./Data/Water Bodies/distToNearestCoastMeters_utm38s.tif')
				# rast <- clusterR(rast, log10Fx)
				# staticPred <- predictAndSumFx(thisPred=staticPred, rast=rast, beta=beta)

			# }

			# # distance to nearest inland water body
			# if ('distToInlandWater' %in% variables) {

				# say('| distToInlandWater', post=0)

				# beta <- betaHats[['distToInlandWater']]

				# rast <- raster('./Data/Water Bodies/distToNearestInlandWaterMeters_utm38s.tif')
				# rast <- clusterR(rast, log10Fx)
				# staticPred <- predictAndSumFx(thisPred=staticPred, rast=rast, beta=beta)

			# }

			# # un/protected
			# if ('as.factor(protectedArea)1' %in% variables) {

				# say('| protectedArea')

				# beta <- betaHats[['as.factor(protectedArea)1']]

				# rast <- raster('./Data/Protected Areas/wdpa_utm38s.tif')
				# rast <- clusterR(rast, calc, args=list(fun=naToZeroFx))

				# staticPred <- predictAndSumFx(thisPred=staticPred, rast=rast, beta=beta)

			# }

		# names(staticPred) <- 'defoLocationStaticPredictors_utm38s'
		# writeRaster(staticPred, './Deforestation Models/Deforestation Location Probability 2015-2080/!defoLocationStaticPredictors_utm38s')

	# } # if needing de novo static predictions

	# ### predict by year using dynamic predictors
	# ############################################


	# say('Predicting using dynamic predictors...')

	# humidForestBufferMask_utm38s <- raster('./Study Region & Masks/UTM 38S 30-m Resolution/humidForestBufferMask_utm38s.tif')

	# for (predictToYear in 2015:2080) {

		# tictoc::tic()
		# say('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++', pre=2)
		# say('Predicting using dynamic predictors for ', predictToYear, ':', post=0)

		# dynPred <- staticPred

		# forestFrom <- if (predictToYear == 2015) {
			# raster('./Data/Forest - Vieilledent et al 2018/forest2014.tif')
		# } else {
			# raster(paste0('./Deforestation Models/Forest 2015-2080 Assuming Pop-Sensitive Loss Amount/forest', predictToYear - 1, '.tif'))
		# }

		# ### distance to forest edge
		# ###########################

		# if ('distToForestEdge' %in% names(betaHats)) {

			# say('dist to forest edge', post=0)

			# if (predictToYear == 2015) {
				# rast <- raster('./Data/Forest - Vieilledent et al 2018/distToForestEdgeMeters2014_utm38s.tif')
			# } else {

				# rast <- fasterRastDistance(forestFrom, metric='euclidean', meters=TRUE, fillNAs=FALSE, grassDir=grassDir)
				# rast <- round(rast)
				# rast <- rast * humidForestBufferMask_utm38s
				# names(rast) <- 'distToForestEdgeMeters2014_utm38s'

				# dirCreate('./Deforestation Models/Distance To Forest Edge 2015-2080')

				# writeRaster(rast, paste0('./Deforestation Models/Distance To Forest Edge 2015-2080/distToForestEdgeMeters', predictToYear, '_utm38s'), datatype='INT4S')

			# }

			# rast <- clusterR(rast, log10Fx)

			# beta <- betaHats[['distToForestEdge']]
			# dynPred <- predictAndSumFx(thisPred=dynPred, rast=rast, beta=beta)

		# }

		# ### distance to deforestation
		# #############################

		# if ('distToDefo' %in% names(betaHats)) {

			# say('| dist to deforestation', post=0)

			# if (predictToYear == 2015) {
				# rast <- raster('./Data/Forest - Vieilledent et al 2018/distToDefoMeters2014_utm38s.tif')
			# } else {

				# # calculate distance to most recent deforestation using prior year's forest and forest from year prior to prior year
				# forestBeforeFrom <- if (predictToYear == 2016) {
					# raster('./Data/Forest - Vieilledent et al 2018/forest2014.tif')
				# } else if (predictToYear > 2016) {
					# raster(paste0('./Deforestation Models/Forest 2015-2080 Assuming Pop-Sensitive Loss Amount/forest', predictToYear - 2, '.tif'))
				# }

				# # identify areas of deforestation
				# forestFromZeros <- clusterR(forestFrom, calc, args=list(fun=naToZeroFx))

				# defoFromVsBeforeFrom <- forestBeforeFrom - forestFromZeros
				# defoFromVsBeforeFrom <- clusterR(defoFromVsBeforeFrom, calc, args=list(fun=nonOnesToNaFx))

				# rast <- fasterRastDistance(defoFromVsBeforeFrom, metric='euclidean', meters=TRUE, fillNAs=TRUE, grassDir=grassDir)

				# rast <- round(rast)
				# rast <- rast * humidForestBufferMask_utm38s
				# names(rast) <- paste0('distToDefoMeters', predictToYear, '_utm38s')

				# dirCreate('./Deforestation Models/Distance To Deforestation 2015-2080')

				# writeRaster(rast, paste0('./Deforestation Models/Distance To Deforestation 2015-2080/distToDefoMeters', predictToYear, '_utm38s'), datatype='INT4S')

			# }

			# rast <- clusterR(rast, log10Fx)

			# beta <- betaHats[['distToDefo']]
			# dynPred <- predictAndSumFx(thisPred=dynPred, rast=rast, beta=beta)

		# }

		# ### forest fragmentation class
		# ##############################

		# if ('as.factor(fragClass)6' %in% names(betaHats)) {

			# say('| fragClass', post=0)

			# if (predictToYear == 2015) {
				# fragClass <- raster('./Data/Forest - Vieilledent et al 2018/forestFragClass2014_utm38s.tif')
			# } else {

				# forestFromZeros <- clusterR(forestFrom, calc, args=list(fun=naToZeroFx))

				# frag <- fasterFragmentation(rast = forestFromZeros, size = 5, pad = TRUE, padValue = NA, calcDensity = TRUE, calcConnect = TRUE, calcClass = TRUE, na.rm = TRUE, undet = 'edge', cores = 4, forceMulti = TRUE)

				# fragClass <- frag[['class']]
				# fragClass <- humidForestBufferMask_utm38s * frag
				# names(fragClass) <- paste0('forestFragClass', predictToYear - 1, '_utm38s')

				# dirCreate('./Deforestation Models/Fragmentation Class 2015-2080')

				# writeRaster(fragClass, paste0('./Deforestation Models/Fragmentation Class 2015-2080/forestFragClass', predictToYear - 1, '_utm38s'), datatype='INT1U')

			# }

			# fiveToFourFunct <- function(x) ifelse(x == 5, 4, x)
			# fragClass <- clusterR(fragClass, calc, args=list(fun=fiveToFourFunct))

			# onlyThisLevel <- function(x, lev) ifelse(x == lev, 1, NA)

			# # predict for each level (note level 1 is assumed in intercept)
			# for (thisLevel in c(2, 3, 4, 6)) {

				# rast <- clusterR(fragClass, fun=onlyThisLevel, args=list(lev=thisLevel))
				# beta <- betaHats[[paste0('as.factor(fragClass)', thisLevel)]]
				# dynPred <- predictAndSumFx(thisPred=dynPred, rast=rast, beta=beta)

			# }

		# }

		# ### predict deforestation location
		# ##################################

			# say('| deforestation location probability', post=0)

			# invLogit <- function(x, ...) phcfM::inv.logit(x)
			# probOfDefo <- clusterR(dynPred, calc, args=list(fun=invLogit))
			# probOfDefo <- probOfDefo * forestFrom

			# names(probOfDefo) <- paste0('defoLocationProb', predictToYear - 1, 'to', predictToYear)

			# dirCreate('./Deforestation Models/Deforestation Location Probability 2015-2080')
			# writeRaster(probOfDefo, paste0('./Deforestation Models/Deforestation Location Probability 2015-2080/defoLocationProb', predictToYear - 1, 'to', predictToYear))

		# ### predict forest cover for this year
		# ######################################

			# say('| forest cover', post=0)

			# # ### calculate new forest cover using assuming loss RATE is constant
			# # ###################################################################

				# # thold <- fasterQuantile(probOfDefo, probs=1 - r, grassDir=grassDir)

				# # tholdFx <- function(x) ifelse(x > thold, 1, 0)
				# # thisDefo <- clusterR(probOfDefo, calc, args=list(fun=tholdFx), export='thold')

				# # newForest <- forestFrom - thisDefo
				# # names(newForest) <- paste0('forest', predictToYear, '_utm38s')

				# # dirCreate('./Deforestation Models/Forest 2015-2080 Assuming Constant Loss Rate')
				# # writeRaster(newForest, paste0('./Deforestation Models/Forest 2015-2080 Assuming Constant Loss Rate/forest', predictToYear), datatype='INT1U')

			# ### calculate new forest cover using assuming loss AMOUNT varies with time
			# ##########################################################################

				# # area of forest to lose
				# forestToLose <- forestToLoseBase <- amounts$amount[amounts$year == predictToYear]

				# # amount of forest remaining in year prior
				# forestLeft <- cellStats(forestFrom, 'sum')

				# if (forestLeft < forestToLose) {
					# newForest <- NA * forestFrom
				# } else {

					# # adjust threshold if amount lost is different from amount desired to be lost (note sure why this happens... seems to be something about the quantile function not choosing exactly the right value

					# ratioLostToDesiredToLose <- 1

					# say('tries:', post=0)
					# tries <- 1
					# while (tries == 1 | (ratioLostToDesiredToLose > 1.01 | ratioLostToDesiredToLose < 0.99)) {

						# say(tries, post=0)
						# tries <- tries + 1

						# forestToLose <- forestToLose / ratioLostToDesiredToLose

						# # get threshold for deciding what's deforested/forested
						# probs <- (forestLeft - forestToLose) / forestLeft
						# thold <- fasterQuantile(probOfDefo, probs=probs, grassDir=grassDir)

						# # apply threshold
						# tholdFx <- function(x) ifelse(x > thold, NA, 1)
						# thisNonDefo <- clusterR(probOfDefo, calc, args=list(fun=tholdFx), export='thold')

						# newForestLeft <- cellStats(thisNonDefo, 'sum')

						# forestActuallyLost <- forestLeft - newForestLeft
						# ratioLostToDesiredToLose <- forestActuallyLost / forestToLoseBase

						# say('forestLeft: ', forestLeft, ' forestToLose: ', forestToLose, ' forestActuallyLost: ', forestActuallyLost, ' ratioLostToDesiredToLose: ', ratioLostToDesiredToLose, ' probs: ', probs, ' thold: ', thold, pre=1)

					# }

				# }

				# newForest <- forestFrom * thisNonDefo
				# names(newForest) <- paste0('forest', predictToYear, '_utm38s')

				# dirCreate('./Deforestation Models/Forest 2015-2080 Assuming Pop-Sensitive Loss Amount')
				# writeRaster(newForest, paste0('./Deforestation Models/Forest 2015-2080 Assuming Pop-Sensitive Loss Amount/forest', predictToYear), datatype='INT1U')

			# rm(dynPred, probOfDefo, forestFrom, newForest, thisNonDefo)
			# gc()

			# say('| ', post=0)
			# tictoc::toc()

	# } # next predict-to year

	# endCluster()

	
# say('####################################################################################################')
# say('### create maps of future forest cover assuming no additional loss in protected areas after 2014 ###')
# say('####################################################################################################')

	# say('Assuming forest cover stays unchanged from 2014 in protected areas and that foregone forest products from these areas are do not displace demand to other areas in Madagascar.')

	# # create mask with current forest cover in PAs
	# pas_utm38s <- raster('./Data/Protected Areas/wdpa_utm38s.tif')
	
	# forest2014 <- raster('./Data/Forest - Vieilledent et al 2018/forest2014.tif')
	# pasWithCurrentForest_utm38s <- forest2014 * pas_utm38s
	
	# # mask future forest with current PA forest cover
	# forest2050 <- raster('./Deforestation Models/Forest 2015-2080 Assuming Pop-Sensitive Loss Amount/forest2050.tif')
	# forest2070 <- raster('./Deforestation Models/Forest 2015-2080 Assuming Pop-Sensitive Loss Amount/forest2070.tif')
	
	# forest2050pa <- stack(forest2050, pasWithCurrentForest_utm38s)
	# forest2070pa <- stack(forest2070, pasWithCurrentForest_utm38s)

	# forest2050pasMasked_utm38s <- max(forest2050pa, na.rm=TRUE)
	# forest2070pasMasked_utm38s <- max(forest2070pa, na.rm=TRUE)
		
	# beginCluster(6)

		# naFx <- function(x) ifelse(x == 0, NA, x)
		# forest2050pasMasked_utm38s <- clusterR(forest2050pasMasked_utm38s, calc, args=list(fun=naFx))
		# forest2070pasMasked_utm38s <- clusterR(forest2070pasMasked_utm38s, calc, args=list(naFx))
		
	# endCluster()
	
	# names(forest2050pasMasked_utm38s) <- 'forest2050_utm38s'
	# names(forest2070pasMasked_utm38s) <- 'forest2070_utm38s'
	
	# dirCreate('./Deforestation Models/Forest 2015-2080 Assuming Pop-Sensitive Loss Amount PAs Have 2014 Cover')
	# writeRaster(forest2050pasMasked_utm38s, './Deforestation Models/Forest 2015-2080 Assuming Pop-Sensitive Loss Amount PAs Have 2014 Cover/forest2050', dtatype='INT1U')
	# writeRaster(forest2070pasMasked_utm38s, './Deforestation Models/Forest 2015-2080 Assuming Pop-Sensitive Loss Amount PAs Have 2014 Cover/forest2070', dtatype='INT1U')

# say('###################################################')
# say('### calculate future forest fragmentation class ###')
# say('###################################################')

	# cores <- 4

	# zeroFunct <- function(x) ifelse(is.na(x), 0, x)
	# humidForestBufferMask_utm38s <- raster('./Study Region & Masks/UTM 38S 30-m Resolution/humidForestBufferMask_utm38s.tif')
	
	# # forest loss anywhere
	# for (year in c(2050, 2070)) {
	
		# say('forest loss anywhere: ', year)
	
		# forest <- raster(paste0('./Deforestation Models/Forest 2015-2080 Assuming Pop-Sensitive Loss Amount/forest', year, '.tif'))
	
		# beginCluster(7)
			# forestZeros <- clusterR(forest, calc, args=list(fun=zeroFunct))
		# endCluster()
		
		# frag <- fasterFragmentation(rast = forestZeros, size = 5, pad = TRUE, padValue = NA, calcDensity = TRUE, calcConnect = TRUE, calcClass = TRUE, na.rm = TRUE, undet = 'perforated', cores = cores, forceMulti = TRUE)

		# frag <- frag[['class']]
		# frag <- frag * humidForestBufferMask_utm38s
		
		# names(frag) <- paste0('forestFragClass', year, '_utm38s')

		# dirCreate('./Deforestation Model/Forest Fragmentation Class 2015-2080 Assuming Pop-Sensitive Loss Amount')
		# writeRaster(frag, paste0('./Deforestation Model/Forest Fragmentation Class 2015-2080 Assuming Pop-Sensitive Loss Amount/forestFragClass', year, '_utm38s'), datatype='INT1U')
		
		# rm(frag); gc()
		
	# } # next year
	
	# # forest loss anywhere but protected areas
	# for (year in c(2050, 2070)) {
	
		# say('forest loss anywhere but protected: ', year)
	
		# forest <- raster(paste0('./Deforestation Models/Forest 2015-2080 Assuming Pop-Sensitive Loss Amount PAs Have 2014 Cover/forest', year, '.tif'))
	
		# beginCluster(7)
			# forestZeros <- clusterR(forest, calc, args=list(fun=zeroFunct))
		# endCluster()
		
		# frag <- fasterFragmentation(rast = forestZeros, size = 5, pad = TRUE, padValue = NA, calcDensity = TRUE, calcConnect = TRUE, calcClass = TRUE, na.rm = TRUE, undet = 'perforated', cores = cores, forceMulti = TRUE)

		# frag <- frag[['class']]
		# frag <- frag * humidForestBufferMask_utm38s
		
		# names(frag) <- paste0('forestFragClass', year, '_utm38s')

		# dirCreate('./Deforestation Model/Forest Fragmentation Class 2015-2080 Assuming Pop-Sensitive Loss Amount PAs Have 2014 Cover')
		# writeRaster(frag, paste0('./Deforestation Model/Forest Fragmentation Class 2015-2080 Assuming Pop-Sensitive Loss Amount PAs Have 2014 Cover/forestFragClass', year, '_utm38s'), datatype='INT1U')
		
		# rm(frag); gc()
		
	# } # next year

# say('#########################################################')
# say('### compile forest cover and fragmentation statistics ###')
# say('#########################################################')

	# forest2014 <- raster('./Data/Forest - Vieilledent et al 2018/forest2014.tif')
	# frag2014 <- raster('./Data/Forest - Vieilledent et al 2018/forestFragClass2014_utm38s.tif')
	
	# forest2050 <- raster('./Deforestation Models/Forest 2015-2080 Assuming Pop-Sensitive Loss Amount/forest2050.tif')
	# frag2050 <- raster('./Deforestation Models/Forest Fragmentation Class 2015-2080 Assuming Pop-Sensitive Loss Amount/forestFragClass2050_utm38s.tif')
	
	# forest2070 <- raster('./Deforestation Models/Forest 2015-2080 Assuming Pop-Sensitive Loss Amount/forest2070.tif')
	# frag2070 <- raster('./Deforestation Models/Forest Fragmentation Class 2015-2080 Assuming Pop-Sensitive Loss Amount/forestFragClass2070_utm38s.tif')
	
	# forest2050pa <- raster('./Deforestation Models/Forest 2015-2080 Assuming Pop-Sensitive Loss Amount PAs Have 2014 Cover/forest2050.tif')
	# frag2050pa <- raster('./Deforestation Models/Forest Fragmentation Class 2015-2080 Assuming Pop-Sensitive Loss Amount PAs Have 2014 Cover/forestFragClass2050_utm38s.tif')
	
	# forest2070pa <- raster('./Deforestation Models/Forest 2015-2080 Assuming Pop-Sensitive Loss Amount PAs Have 2014 Cover/forest2070.tif')
	# frag2070pa <- raster('./Deforestation Models/Forest Fragmentation Class 2015-2080 Assuming Pop-Sensitive Loss Amount PAs Have 2014 Cover/forestFragClass2070_utm38s.tif')
	
	# out <- data.frame()
	
	# for (year in c(2014, 2050, 2070)) {
	
		# say(year)
	
		# forest <- get(paste0('forest', year))
		# forest <- forest > 0
		# forest <- cellStats(forest, 'sum')
		# forest <- (30^2 * forest) / 1000^2
	
		# interior <- get(paste0('frag', year))
		# interior <- interior == 6
		# interior <- cellStats(interior, 'sum')
		# interior <- (30^2 * interior) / 1000^2
		
		# out <- rbind(
			# out,
			# data.frame(
				# year = year,
				# forest_km2 = forest,
				# interior_km2 = interior,
				# pa = FALSE
			# )
		# )
		
		# # defo in protected areas prohibited
		# if (year > 2014) {
			
			# say(year, ' pa')
			
			# forest <- get(paste0('forest', year, 'pa'))
			# forest <- forest > 0
			# forest <- cellStats(forest, 'sum')
			# forest <- (30^2 * forest) / 1000^2
		
			# interior <- get(paste0('frag', year, 'pa'))
			# interior <- interior == 6
			# interior <- cellStats(interior, 'sum')
			# interior <- (30^2 * interior) / 1000^2
			
			# out <- rbind(
				# out,
				# data.frame(
					# year = year,
					# forest_km2 = forest,
					# interior_km2 = interior,
					# pa = TRUE
				# )
			# )
		
		# } # PA
		
	# }

	# dirCreate('./Figures & Tables/Forest - Change Statistics')
	# write.csv(out, './Figures & Tables/Forest - Change Statistics/Forest Area and Interior Area by Scenario.csv', row.names=FALSE)


# say('###########################################')
# say('### create display maps of forest cover ###')
# say('###########################################')

	# # generalization
	# outDir <- './Figures & Tables/Forest - Cover/'

	# humidForestBufferMask_utm38s <- raster('./Study Region & Masks/UTM 38S 30-m Resolution/humidForestBufferMask_utm38s.tif')

	# load('./Study Region & Masks/UTM 38S 30-m Resolution/Eastern Humid Forest Polygon.RData')
	# load('./Study Region & Masks/UTM 38S 30-m Resolution/Madagascar from GADM 3.6.RData')
	
	# # pas <- shapefile('./Data/Protected Areas/WDPA_Sept2018_MDG-shapefile-polygons')
	# # pas <- sp::spTransform(pas, CRS(madEaProj))
	# # pas <- crop(pas, madagascar_utm38s)
	
	# load('./Data/Protected Areas/WDPA_Sept2018_MDG-shapefile-polygons-onlyTerrestrial.RData')
	
	# cols <- 'forestgreen'
	# names(cols) <- 'Forest'
	
	# pars <- par()

	# for (year in c(2014, 2050, 2070)) {
	
		# defos <- if (year <= 2014) { 'anywhere' } else {c('anywhere', 'notPAs') }
	
		# for (defo in defos) {
		
			# say(year, ' ', defo)
			
			# if (year <= 2014) {
				# x <- raster(paste0('./Data/Forest - Vieilledent et al 2018/forest', year, '.tif'))
				# pngName <- paste0('Forest Cover for ', year)
			# } else if (defo == 'anywhere') {
				# x <- raster(paste0('./Deforestation Models/Forest 2015-2080 Assuming Pop-Sensitive Loss Amount/forest', year, '.tif'))
				# pngName <- paste0('Forest Cover for ', year, ' Assuming Relaxed Protection of Forest')
			# } else if (defo == 'notPAs') {
				# x <- raster(paste0('./Deforestation Models/Forest 2015-2080 Assuming Pop-Sensitive Loss Amount PAs Have 2014 Cover/forest', year, '.tif'))
				# pngName <- paste0('Forest Cover for ', year, ' for Assuming PAs have 2014 Cover')
			# }
			
			# dirCreate(outDir)
			
			# # all of Madagascar
			# png(paste0(outDir, '/', pngName, '.png'), width=400, height=1200, res=200)

				# par(mar=0.1 * rep(1, 4), oma=0.1 * rep(1, 4))

				# # eastern moist forest
				# plot(humidForest_utm38s, lwd=0.5)
				# plot(x, col=cols, legend=FALSE, add=TRUE)
				# plot(pas, border='blue', add=TRUE)
				# plot(madagascar_utm38s, xpd=NA, lwd=0.5, add=TRUE)

				# # insets
				# plot(madFocus1, lwd=2, add=TRUE)
				# plot(madFocus2, lwd=2, add=TRUE)
				# plot(madFocus3, lwd=2, add=TRUE)

				# legend('bottomright', inset=c(0, 0.05), legend=c(names(cols), 'Protected area'), fill=c(cols, NA), border=c(rep('black', length(cols)), 'blue'), cex=0.7, bty='n')
				
			# dev.off()

		# } # next defo location
		
	# } # next year
	
	# par(pars)
	
# say('#############################################')
# say('### compile GIS files for faster plotting ###')
# say('#############################################')
	
	# # humid forest shapefile
	# load('./Study Region & Masks/UTM 38S 30-m Resolution/Eastern Humid Forest Polygon.RData')

	# ### MOBOT PAs
	# #############
	
	# mobotPaFiles <- listFiles('./Data/Protected Areas/MOBOT/shapefiles', pattern='.shp')
	# mobotPaFiles <- mobotPaFiles[!grepl(mobotPaFiles, pattern='.xml')]
	# mobotPas <- list()
	# for (mobotPaFile in mobotPaFiles) {
		# thisPa <- shapefile(mobotPaFile)
		# if (is.na(projection(thisPa))) projection(thisPa) <- getCRS('wgs84')
		# thisPa <- sp::spTransform(thisPa, CRS(madEaProj))
		# inHumidForest <- over(thisPa, humidForest_utm38s)
		# if (!all(is.na(inHumidForest))) {
			# thisPa <- as(thisPa, 'SpatialPolygons')
			# projection(thisPa) <- madEaProj
			# mobotPas[[length(mobotPas) + 1]] <- thisPa
			# name <- basename(mobotPaFile)
			# name <- gsub(pattern='.shp', name, replacement='')
			# names(mobotPas)[length(mobotPas)] <- name
		# }
	# }
	
	# # combine MOBOT PAs
	# mobotPasCombine <- list(mobotPas, makeUniqueIDs = T) %>%
		# flatten() %>%
		# do.call(rbind, .)
		
	# save(mobotPas, file='./Data/Protected Areas/MOBOT/MOBOT PAs in Eastern Humid Forest as List.RData')
	# save(mobotPasCombine, file='./Data/Protected Areas/MOBOT/MOBOT PAs in Eastern Humid Forest as SpatialPolygon.RData')
	
	# # shapefile(mobotPasCombine, 'C:/ecology/!Scratch/mobotPas') # for rasterization in next step
	
	# # ### raster mask for PAs including MOBOT's
	# # pas_utm38s <- raster('./Data/Protected Areas/wdpa_utm38s.tif')
	# # mobotPas <- raster('**************') # generated in QGIS

	# # p <- stack(pas_utm38s, r)
	# # beginCluster(6)
		# # fx <- function(x) ifelse(any(!is.na(x)), 1, NA)
		# # x0InPas <- clusterR(p, calc, args=list(fun=fx))
	# # endCluster()
	
	# # pasWpdaMobot_utm838s <- x0InPas
	# # pasWpdaMobot_utm838s <- round(pasWpdaMobot_utm838s)
	# # pasWpdaMobot_utm838s <- setMinMax(pasWpdaMobot_utm838s)
	# # names(pasWpdaMobot_utm838s) <- 'pasWpdaMobot_utm838s'
	# # writeRaster(pasWpdaMobot_utm838s, './Data/Protected Areas/wdpaMobot_utm38s', datatype='INT1U')
	
	# ### forest in PAS in 2014
	# x0 <- raster('./Data/Forest - Vieilledent et al 2018/forest2014.tif')
	
	# forest2014InPas <- x0 * pasWpdaMobot_utm838s
	# forest2014InPas <- round(forest2014InPas)
	# forest2014InPas <- setMinMax(forest2014InPas)
	# projection(forest2014InPas) <- madEaProj
	# writeRaster(forest2014InPas, './Data/Forest - Vieilledent et al 2018/forest2014InPAsWpdaMobot', datatype='INT1U')
	
	# ### eastern humid forest mask (portion OUTSIDE the buffer)
	
	# # humid forest shapefile
	# hs_utm38s <- raster('./Data/Topography - GMTED2010/hillshadeGmted2010_utm38s.tif')
	# ext <- extent(hs_utm38s)
	# ext <- as(ext, 'SpatialPolygons')
	# projection(ext) <- madEaProj

	# load('./Study Region & Masks/UTM 38S 30-m Resolution/Eastern Humid Forest Polygon.RData')
	
	# notHumidForest_utm38s <- gDifference(ext, humidForest_utm38s)
	# save(notHumidForest_utm38s, file='./Study Region & Masks/UTM 38S 30-m Resolution/NOT Eastern Humid Forest Polygon.RData')

# say('################################################################')
# say('### create display maps of just Madagascar for presentations ###')
# say('################################################################')

	# # generalization
	# outDir <- './Figures & Tables/Forest - Cover for Presentations/'

	# load('./Study Region & Masks/UTM 38S 30-m Resolution/Madagascar from GADM 3.6.RData')
	
	# # humid forest shape
	# load('./Study Region & Masks/UTM 38S 30-m Resolution/NOT Eastern Humid Forest Polygon.RData')

	# ### PAs
	# load('./Data/Protected Areas/WDPA_Sept2018_MDG-shapefile-polygons-onlyTerrestrial.RData')

	# ### MOBOT PAs
	# load('./Data/Protected Areas/MOBOT/MOBOT PAs in Eastern Humid Forest as SpatialPolygon.RData')

	# ### hillshade colors
	# hs <- raster('./Data/Topography - GMTED2010/hillshadeGmted2010_utm38s.tif')
	# hsCols <- paste0('gray', 0:100)
	# hsCols <- alpha(hsCols, 0.8)
	
	# x0col <- alpha('red', 0.6)
	# x1col <- alpha('green', 0.6)
	
	# paBorder <- 'blue'
	# mobotPaBorder <- 'cyan'
	
	# # all of Madagascar
	# png(paste0(outDir, '/Madagascar with Eastern Humid Forest Highlighted.png'), width=1900, height=2400, res=600)

		# par(mar=rep(0, 4), oma=rep(0, 4), bg='black')

		# plot(madagascar_utm38s, border=NA, ann=FALSE)
		# plot(hs, col=hsCols, legend=FALSE, ann=FALSE, add=TRUE)
		# plot(notHumidForest_utm38s, col=alpha('black', 0.4), border=NA, add=TRUE, ann=FALSE)
		# # plot(pas, border=paBorder, lwd=0.5, add=TRUE, ann=FALSE)
		# # plot(madagascar_utm38s, xpd=NA, lwd=0.8, add=TRUE, border=NA, ann=FALSE)
		# plot(mobotPasCombine, add=TRUE, border=mobotPaBorder, lwd=0.4)

	# dev.off()

	# # all of Madagascar
	# png(paste0(outDir, '/Madagascar with Eastern Humid Forest Highlighted and PAs.png'), width=1900, height=2400, res=600)

		# par(mar=rep(0, 4), oma=rep(0, 4), bg='black')

		# plot(madagascar_utm38s, border=NA, ann=FALSE)
		# plot(hs, col=hsCols, legend=FALSE, ann=FALSE, add=TRUE)
		# plot(notHumidForest_utm38s, col=alpha('black', 0.4), border=NA, add=TRUE, ann=FALSE)
		# plot(pas, border=paBorder, lwd=0.5, add=TRUE, ann=FALSE)
		# # plot(madagascar_utm38s, xpd=NA, lwd=0.8, add=TRUE, border=NA, ann=FALSE)
		# plot(mobotPasCombine, add=TRUE, border=mobotPaBorder, lwd=0.4)

	# dev.off()

# say('##################################################################')
# say('### create display maps of PAST forest cover for presentations ###')
# say('##################################################################')

	# # generalization
	# outDir <- './Figures & Tables/Forest - Cover for Presentations/'

	# load('./Study Region & Masks/UTM 38S 30-m Resolution/Madagascar from GADM 3.6.RData')
	
	# # humid forest shape
	# load('./Study Region & Masks/UTM 38S 30-m Resolution/NOT Eastern Humid Forest Polygon.RData')

	# ### PAs
	# load('./Data/Protected Areas/WDPA_Sept2018_MDG-shapefile-polygons-onlyTerrestrial.RData')

	# ### MOBOT PAs
	# load('./Data/Protected Areas/MOBOT/MOBOT PAs in Eastern Humid Forest as SpatialPolygon.RData')

	# ### hillshade colors
	# hs <- raster('./Data/Topography - GMTED2010/hillshadeGmted2010_utm38s.tif')
	# hsCols <- paste0('gray', 0:100)
	# hsCols <- alpha(hsCols, 0.8)
	
	# x0col <- alpha('red', 0.6)
	# x1col <- alpha('green', 0.6)
	
	# paBorder <- 'blue'
	# mobotPaBorder <- 'cyan'
	
	# # for (year in c(2000, 2005, 2010, 2014)) {
	# # for (year in c(2000, 2014)) {
	# for (year in c(2014)) {

		# say(year)

		# if (year != 2000) x0 <- raster(paste0('./Data/Forest - Vieilledent et al 2018/deforest', year, '_utm38s.tif'))
		# x1 <- raster(paste0('./Data/Forest - Vieilledent et al 2018/forest', year, '.tif'))

		# pngName <- paste0('Forest Cover for XX', year)
		# dirCreate(outDir)
		
		# # all of Madagascar
		# png(paste0(outDir, '/', pngName, '.png'), width=1900, height=2400, res=600)

			# par(mar=rep(0, 4), oma=rep(0, 4), bg='black')

			# plot(madagascar_utm38s, border=NA, ann=FALSE)
			# plot(hs, col=hsCols, legend=FALSE, ann=FALSE, add=TRUE)
			# if (year != 2000) plot(x0, col=x0col, legend=FALSE, add=TRUE)
			# plot(x1, col=x1col, legend=FALSE, add=TRUE)
			# plot(notHumidForest_utm38s, col=alpha('black', 0.4), border=NA, add=TRUE, ann=FALSE)
			# # plot(pas, border=paBorder, lwd=0.5, add=TRUE, ann=FALSE)
			# # plot(madagascar_utm38s, xpd=NA, lwd=0.8, add=TRUE, border=NA, ann=FALSE)
			# plot(mobotPasCombine, add=TRUE, border=mobotPaBorder, lwd=0.4)

			# text(470000, 8573849, labels=year, cex=1.4, xpd=NA, pos=4, col='white')

		# dev.off()

	# } # next year

# say('####################################################################')
# say('### create display maps of FUTURE forest cover for presentations ###')
# say('####################################################################')

	# # generalization
	# outDir <- './Figures & Tables/Forest - Cover for Presentations/'

	# # humid forest shape
	# load('./Study Region & Masks/UTM 38S 30-m Resolution/NOT Eastern Humid Forest Polygon.RData')

	# # Madagascar shapefile
	# load('./Study Region & Masks/UTM 38S 30-m Resolution/Madagascar from GADM 3.6.RData')
	
	# ### PAs
	# load('./Data/Protected Areas/WDPA_Sept2018_MDG-shapefile-polygons-onlyTerrestrial.RData')
	# load('./Data/Protected Areas/MOBOT/MOBOT PAs in Eastern Humid Forest as List.RData')
	# load('./Data/Protected Areas/MOBOT/MOBOT PAs in Eastern Humid Forest as SpatialPolygon.RData')

	# ### 2014 forest cover raster in PAs
	# forest2014InPAsWpdaMobot <- raster('./Data/Forest - Vieilledent et al 2018/forest2014InPAsWpdaMobot.tif')

	# ### hillshade
	# hs <- raster('./Data/Topography - GMTED2010/hillshadeGmted2010_utm38s.tif')
	# hsCols <- paste0('gray', 0:100)
	# hsCols <- alpha(hsCols, 0.8)
	
	# x0col <- alpha('red', 0.6)
	# x1col <- alpha('green', 0.6)
	
	# paBorder <- 'blue'
	# mobotPaBorder <- 'cyan'
	
	# for (year in c(2015:2080)) {
	# # for (year in c(2015)) {

		# for (protection in c('anywhere', 'notPAs')) {
		
			# say(year, ' ', protection)

			# x1 <- raster(paste0('./Deforestation Models/Forest 2015-2080 Assuming Pop-Sensitive Loss Amount/forest', year, '.tif'))
			
			# if (protection == 'anywhere') {
				
				# pngName <- paste0('Forest Cover for ', year, ' - Deforestation Anywhere')
				# thisOutDir <- paste0(outDir, '/Relaxed Protection')

			# } else if (protection == 'notPAs') {
			
				# pngName <- paste0('Forest Cover for ', year, ' - PAs have 2014 Cover')
				# thisOutDir <- paste0(outDir, '/Strict Protection')

				# # add forest in PAs
				# x1pas <- stack(x1, forest2014InPAsWpdaMobot)
				# x1pas <- max(x1pas, na.rm=TRUE)
					
				# beginCluster(6)

					# fx <- function(x) ifelse(x == 1, 1, NA)
					# x1 <- clusterR(x1pas, calc, args=list(fun=fx))
					
				# endCluster()
				
			# }

			# dirCreate(thisOutDir)
			
			# # all of Madagascar
			# png(paste0(thisOutDir, '/', pngName, '.png'), width=1900, height=2400, res=600)

				# par(mar=rep(0, 4), oma=rep(0, 4), bg='black')

				# # eastern moist forest
				# plot(madagascar_utm38s, border=NA, ann=FALSE)
				# plot(hs, col=hsCols, legend=FALSE, ann=FALSE, add=TRUE)
				# # plot(humidForest_utm38s, lwd=0.2, ann=FALSE, add=TRUE)
				# # plot(x0mask, col=x0col, legend=FALSE, add=TRUE)
				# plot(x1, col=x1col, legend=FALSE, add=TRUE)
				# plot(notHumidForest_utm38s, col=alpha('black', 0.4), border=NA, add=TRUE, ann=FALSE)
				# if (protection == 'notPAs') plot(pas, border=paBorder, lwd=0.4, add=TRUE, ann=FALSE)
				# plot(mobotPasCombine, add=TRUE, border=mobotPaBorder, lwd=0.4)

				# text(470000, 8573849, labels=year, cex=1.4, xpd=NA, pos=4, col='white')

			# dev.off()

		# } # next forest protection

	# } # next year

# say('###############################################################################')
# say('### create display maps of forest cover in MOBOT reserves for presentations ###')
# say('###############################################################################')

	# # humid forest shape
	# load('./Study Region & Masks/UTM 38S 30-m Resolution/NOT Eastern Humid Forest Polygon.RData')

	# ### MOBOT PAs
	# load('./Data/Protected Areas/MOBOT/MOBOT PAs in Eastern Humid Forest as List.RData')
	# load('./Data/Protected Areas/MOBOT/MOBOT PAs in Eastern Humid Forest as SpatialPolygon.RData')

	# ### hillshade
	# hs <- raster('./Data/Topography - GMTED2010/hillshadeGmted2010_utm38s.tif')
	# hsCols <- paste0('gray', 0:100)
	# hsCols <- alpha(hsCols, 0.8)
	
	# x0col <- alpha('red', 0.6)
	# x1col <- alpha('green', 0.6)

	# # elevation
	# elev <- raster('./Data/Topography - GMTED2010/elevationGmted2010_utm38s.tif')
	# elevMinMax <- c(minValue(elev), maxValue(elev))
	# elevRange <- diff(elevMinMax)
	
	# # 2000 forest cover raster
	# x0 <- raster(paste0('./Data/Forest - Vieilledent et al 2018/forest2000.tif'))

	# # 2014 forest cover raster
	# x1 <- raster(paste0('./Data/Forest - Vieilledent et al 2018/forest2014.tif'))

	# beginCluster(4)
	# for (targetPa in names(mobotPas)) {
		
		# say(targetPa)
		
		# thisPa <- mobotPas[[targetPa]]
		# thisPa <- sp::spTransform(thisPa, CRS(projection(x0)))
		# ext <- extent(thisPa)
		# bufferSize <- 1.2 * mean(c(ext@xmax - ext@xmin, ext@ymax - ext@ymin))
		# thisExt <- gBuffer(thisPa, width=bufferSize)
		
		# elevCrop <- crop(elev, thisExt)
		# hsCrop <- crop(hs, thisExt)
		# x0crop <- crop(x0, thisExt)
		# x1crop <- crop(x1, thisExt)

		# # mask to just portions lost
		# x01crop <- x0crop + x1crop

		# # processed in cluster
			# naFx <- function(x) ifelse(is.na(x), 1, NA)
			# x0mask <- clusterR(x01crop, calc, args=list(fun=naFx))
		# # end process in cluster
			
		# x0crop <- x0crop * x0mask
		
		# png(paste0('./Figures & Tables/Forest - Cover for Presentations/MOBOT PAs - ', targetPa, ' 2000-2014.png'), width=1600, height=1600, res=800)
		
			# par(mar=rep(0, 4), oma=rep(0, 4), bg='black')
			# plot(thisExt, border=NA, ann=FALSE)
			
			# plot(hsCrop, col=hsCols, legend=FALSE, ann=FALSE, add=TRUE)
			# plot(notHumidForest_utm38s, col=alpha('black', 0.4), border=NA, add=TRUE, ann=FALSE)
			# plot(x0crop, col=x0col, legend=FALSE, add=TRUE)
			# plot(x1crop, col=x1col, legend=FALSE, add=TRUE)

			# plot(mobotPasCombine, add=TRUE, border='cyan', lwd=0.8)
			
		# dev.off()
		
	# } # next PA
	# endCluster()

# say('####################################################################################################')
# say('### create maps of future forest cover assuming no additional loss in protected areas after 2014 ###')
# say('####################################################################################################')

	# say('Assuming forest cover stays unchanged from 2014 in protected areas and that foregone forest products from these areas are do not displace demand to other areas in Madagascar.')

	# # create mask with current forest cover in PAs
	# pas_utm38s <- raster('./Data/Protected Areas/wdpa_utm38s.tif')
	
	# forest2014 <- raster('./Data/Forest - Vieilledent et al 2018/forest2014.tif')
	# pasWithCurrentForest_utm38s <- forest2014 * pas_utm38s
	
	# # mask future forest with current PA forest cover
	# forest2050 <- raster('./Deforestation Models/Forest 2015-2080 Assuming Pop-Sensitive Loss Amount/forest2050.tif')
	# forest2070 <- raster('./Deforestation Models/Forest 2015-2080 Assuming Pop-Sensitive Loss Amount/forest2070.tif')
	
	# forest2050pa <- stack(forest2050, pasWithCurrentForest_utm38s)
	# forest2070pa <- stack(forest2070, pasWithCurrentForest_utm38s)

	# forest2050pasMasked_utm38s <- max(forest2050pa, na.rm=TRUE)
	# forest2070pasMasked_utm38s <- max(forest2070pa, na.rm=TRUE)
		
	# beginCluster(6)

		# naFx <- function(x) ifelse(x == 0, NA, x)
		# forest2050pasMasked_utm38s <- clusterR(forest2050pasMasked_utm38s, calc, args=list(fun=naFx))
		# forest2070pasMasked_utm38s <- clusterR(forest2070pasMasked_utm38s, calc, args=list(naFx))
		
	# endCluster()
	
	# names(forest2050pasMasked_utm38s) <- 'forest2050_utm38s'
	# names(forest2070pasMasked_utm38s) <- 'forest2070_utm38s'
	
	# dirCreate('./Deforestation Models/Forest 2015-2080 Assuming Pop-Sensitive Loss Amount PAs Have 2014 Cover')
	# writeRaster(forest2050pasMasked_utm38s, './Deforestation Models/Forest 2015-2080 Assuming Pop-Sensitive Loss Amount PAs Have 2014 Cover/forest2050', dtatype='INT1U')
	# writeRaster(forest2070pasMasked_utm38s, './Deforestation Models/Forest 2015-2080 Assuming Pop-Sensitive Loss Amount PAs Have 2014 Cover/forest2070', dtatype='INT1U')
	
	
# say('######################################################')
# say('### create very plain display maps of forest cover ###')
# say('######################################################')

	# # generalization
	# outDir <- './Figures & Tables/Forest - Cover - Very Plain Images/'

	# humidForestBufferMask_utm38s <- raster('./Study Region & Masks/UTM 38S 30-m Resolution/humidForestBufferMask_utm38s.tif')

	# load('./Study Region & Masks/UTM 38S 30-m Resolution/Eastern Humid Forest Polygon.RData')
	# load('./Study Region & Masks/UTM 38S 30-m Resolution/Madagascar from GADM 3.6.RData')
	
	# cols <- 'forestgreen'
	# names(cols) <- 'Forest'
	
	# pars <- par()

	# for (year in c(2014, 2050, 2070)) {
	
		# defos <- if (year <= 2014) { 'anywhere' } else {c('anywhere', 'notPAs') }
	
		# for (defo in defos) {
		
			# say(year, ' ', defo)
			
			# if (year <= 2014) {
				# x <- raster(paste0('./Data/Forest - Vieilledent et al 2018/forest', year, '.tif'))
				# pngName <- paste0('Forest Cover for ', year)
			# } else if (defo == 'anywhere') {
				# x <- raster(paste0('./Deforestation Models/Forest 2015-2080 Assuming Pop-Sensitive Loss Amount/forest', year, '.tif'))
				# pngName <- paste0('Forest Cover for ', year, ' Assuming Relaxed Protection of Forest')
			# } else if (defo == 'notPAs') {
				# x <- raster(paste0('./Deforestation Models/Forest 2015-2080 Assuming Pop-Sensitive Loss Amount PAs Have 2014 Cover/forest', year, '.tif'))
				# pngName <- paste0('Forest Cover for ', year, ' for Assuming PAs have 2014 Cover')
			# }
			
			# dirCreate(outDir)
			
			# # all of Madagascar
			# png(paste0(outDir, '/', pngName, '.png'), width=400, height=1200, res=200)

				# par(mar=0.1 * rep(1, 4), oma=0.1 * rep(1, 4), bg='black', fg='white')

				# # eastern moist forest
				# plot(humidForest_utm38s, lwd=0.5, border=NA)
				# plot(x, col=cols, legend=FALSE, add=TRUE)
				# plot(madagascar_utm38s, border='white', xpd=NA, lwd=0.5, add=TRUE)

				# # insets
				# # plot(madFocus1, lwd=2, add=TRUE)
				# # plot(madFocus2, lwd=2, add=TRUE)
				# plot(madFocus3, lwd=2, add=TRUE)

			# dev.off()

		# } # next defo location
		
	# } # next year
	
	# par(pars)
	
# say('##############################################')
# say('### create 3D display maps of forest cover ###')
# say('##############################################')

	# ### definitions
	# ###############

	# # futYear <- '2050'
	# futYear <- '2070'

	# z <- 7.5

	# # ENM rasters
	# sq <- raster(paste0('./Data/Forest - Vieilledent et al 2018/forest2014.tif'))
	# relaxedProt <- raster(paste0('./Deforestation Models/Forest 2015-2080 Assuming Pop-Sensitive Loss Amount/forest', futYear, '.tif'))
	# strictProt <- raster(paste0('./Deforestation Models/Forest 2015-2080 Assuming Pop-Sensitive Loss Amount PAs Have 2014 Cover/forest', futYear, '.tif'))
	
	# panels <- c('sq', 'strictProt', 'relaxedProt')

	# # elevation
	# elev <- raster('./Data/Topography - GMTED2010/elevationGmted2010_utm38s.tif')
	# projection(elev) <- madEaProj
	
	# ### by inset
	# ############
	
	# for (inset in 1:3) {

		# say(inset, level=2)
	
		# focus <- get(paste0('madFocus', inset))
		# elevCrop <- crop(elev, focus)
		
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

		# ### by forest scenario
		# ######################
		
		# for (panel in panels) {
		
			# say('inset ', inset, ' ', panel)
		
			# forest <- get(panel)
			# forest <- crop(forest, focus)
			
			# beginCluster(4)
				# forest <- clusterR(forest, calc, args=list(fun=naTo0))
			# endCluster()
			
			# # convert forest to 3 layers RGB
			# r <- matrix(
				# raster::extract(forest, raster::extent(forest), method = 'simple'),
				# nrow = ncol(forest),
				# ncol = nrow(forest)
			# )
			
			# r <- 0 * forest
			# g <- forest
			# b <- 0 * r
			# alpha <- g

			# r <- as.matrix(r)
			# g <- as.matrix(g)
			# b <- as.matrix(b)
			# alpha <- as.matrix(alpha)

			# if (anyNA(r)) r[is.na(r)] <- 0
			# if (anyNA(g)) g[is.na(g)] <- 0
			# if (anyNA(b)) b[is.na(b)] <- 0
			# if (anyNA(alpha)) alpha[is.na(alpha)] <- 0

			# alpha <- 0.5 * alpha

			# rgb <- array(c(r, g, b, alpha), c(nrow(elevCrop), ncol(elevCrop), 4))

			# hsForest <- add_overlay(hs, rgb, 4)

			# ## add PA boundaries
			# hsForestPas <- hsForest
			# for (i in 1:nrow(paBorders)) {
				# for (j in 1:ncol(paBorders)) {
					# if (paBorders[i, j] == 1) {
						# hsForestPas[i, j, 1] <- 0
						# hsForestPas[i, j, 2] <- 0
						# hsForestPas[i, j, 3] <- 1
					# }
				# }
			# }

			# ### plot!
			# solidDepth <- min(elevMat, na.rm=TRUE)
			# plot_3d(hsForestPas, elevMat, zscale=z, theta=45, phi=30, water=TRUE, zoom=0.8, fov=60, soliddepth=solidDepth, windowsize=c(1600, 1200), shadow=FALSE, background='white')
			
			# outDir <- paste0('./Figures & Tables/Forest - Cover 3D/Panels for ', futYear)
			# dirCreate(outDir)
			
			# rgl::rgl.snapshot(paste0(outDir, '/Inset ', inset, ' Panel ', panel, '.png'), fmt='png', top=FALSE)
			# rgl::rgl.close()

		# } # next panel/ENM
		
	# } # next inset
	
# say('########################################################')
# say('### create display maps of deforestation probability ###')
# say('########################################################')

	# # masks and country borders
	# load('./Study Region & Masks/UTM 38S 30-m Resolution/Eastern Humid Forest Polygon.RData')
	# load('./Study Region & Masks/UTM 38S 30-m Resolution/Madagascar from GADM 3.6.RData')

	# # # pas <- shapefile('./Data/Protected Areas/WDPA_Sept2018_MDG-shapefile-polygons')
	# # # pas <- sp::spTransform(pas, CRS(madEaProj))
	# # # pas <- crop(pas, madagascar_utm38s)
	
	# load('./Data/Protected Areas/WDPA_Sept2018_MDG-shapefile-polygons-onlyTerrestrial.RData')

	# # graphic parameters
	# cols <- colorRampPalette(colors=c('forestgreen', 'yellow', 'orange', 'red')) # cold to hot
	# cols <- cols(10)

	# pars <- par()

	# ### Madagascar
	# ##############
		
		# for (year in c(2015, 2050, 2070)) {
		
			# say(year)
			# outDir <- './Figures & Tables/Forest - Probability of Deforestation'
			# dirCreate(outDir)
			
			# pngName <- paste0('Probability of Deforestation from ', year - 1 , ' to ', year)

			# # raster with data of interest
			# x <- raster(paste0('./Deforestation Models/Deforestation Location Probability 2015-2080/defoLocationProb', year - 1, 'to', year, '.tif'))

			# # all of Madagascar
			# png(paste0(outDir, '/', pngName, '.png'), width=400, height=1200, res=300)

				# par(mar=0.1 * rep(1, 4), oma=0.1 * rep(1, 4))

				# # eastern moist forest
				# plot(humidForest_utm38s, lwd=0.5)
				# plot(x, col=cols, legend=FALSE, add=TRUE)
				# plot(pas, border='blue', lwd=0.5, add=TRUE)
				# plot(madagascar_utm38s, xpd=NA, lwd=0.5, add=TRUE)

				# # insets
				# plot(madFocus1, lwd=1.2, add=TRUE)
				# plot(madFocus2, lwd=1.2, add=TRUE)
				# plot(madFocus3, lwd=1.2, add=TRUE)

				# labs <- c('0', '0.5', '1')
				# paSwatch <- list(swatchAdjY=c(0, 0.07), col=alpha('blue', 0.3), border='blue', labels='Protected')
				# legendGrad('bottomright', inset=c(0.3, 0.05), width=0.17, height=0.2, labels=labs, cex=0.6, labAdj=0.5, adjX=c(0.2, 0.5), adjY=c(0.19, 0.72), col=cols, title='Likelihood', boxBorder=NA, boxBg=NA, xpd=NA, swatches=list(paSwatch))
				
			# dev.off()
		
		# } # next year
		
		# par(pars)

	# ### insets (one plot for all three insets and all three periods)
	# ################################################################

		# humidForestBufferMask_utm38s <- raster('./Study Region & Masks/UTM 38S 30-m Resolution/humidForestBufferMask_utm38s.tif')
	
		# # list of maps
		# maps <- c(NA, 'ul', 'um', 'ur', 'ml', 'mm', 'mr', 'll', 'lm', 'lr', NA, 'leftCol', 'middleCol', 'rightCol', 'topRow', 'middleRow', 'bottomRow')

		# outDir <- './Figures & Tables/Forest - Probability of Deforestation'
		# dirCreate(outDir)

		# pngName <- paste0('Probability of Deforestation')
		# png(paste0(outDir, '/', pngName, ' Insets.png'), width=3600, height=3600, res=450)
		
			# pars <- par(mfrow=c(4, 4), oma=rep(0, 4), mar=c(0, 0, 0, 0))

			# # by EACH RASTER
			# for (map in maps) {
				
				# say('subplot ', map)
				
				# lab <- if (is.na(map)) {
					# NA
				# } else if (map == 'topRow') {
					# 'Northern Inset'
				# } else if (map == 'middleRow') {
					# paste0('Middle Inset')
				# } else if (map == 'bottomRow') {
					# paste0('Southern Inset')
				# } else if (map == 'leftCol') {
					# paste0('2014-2015')
				# } else if (map == 'middleCol') {
					# paste0('2049-2050')
				# } else if (map == 'rightCol') {
					# paste0('2069-2070')
				# }
				
				# rowCol <- if (is.na(map)) { c(1, 1) } else
					# if (map == 'topRow') { c(2, 1)} else
					# if (map == 'middleRow') { c(3, 1)} else
					# if (map == 'bottomRow') { c(4, 1)} else
					# if (map == 'leftCol') { c(1, 2)} else
					# if (map == 'middleCol') { c(1, 3)} else
					# if (map == 'rightCol') { c(1, 4)} else
					# if (map == 'ul') { c(1, 1) + 1 } else
					# if (map == 'um') { c(1, 2) + 1 } else
					# if (map == 'ur') { c(1, 3) + 1 } else
					# if (map == 'ml') { c(2, 1) + 1 } else
					# if (map == 'mm') { c(2, 2) + 1 } else
					# if (map == 'mr') { c(2, 3) + 1 } else
					# if (map == 'll') { c(3, 1) + 1 } else
					# if (map == 'lm') { c(3, 2) + 1 } else
					# if (map == 'lr') { c(3, 3) + 1 }
				
				# letter <- if (is.na(map)) { NA } else
					# if (map == 'ul') { 'a' } else
					# if (map == 'um') { 'b' } else
					# if (map == 'ur') { 'c' } else
					# if (map == 'ml') { 'd' } else
					# if (map == 'mm') { 'e' } else
					# if (map == 'mr') { 'f' } else
					# if (map == 'll') { 'g' } else
					# if (map == 'lm') { 'h' } else
					# if (map == 'lr') { 'i' }
				
				# inset <- if (is.na(map)) { NA } else
					# if (map == 'ul') { 1 } else
					# if (map == 'um') { 1 } else
					# if (map == 'ur') { 1 } else
					# if (map == 'ml') { 2 } else
					# if (map == 'mm') { 2 } else
					# if (map == 'mr') { 2 } else
					# if (map == 'll') { 3 } else
					# if (map == 'lm') { 3 } else
					# if (map == 'lr') { 3 }
				
				# par(mfg=rowCol)
				# plot.new()
				
				# # top left corner
				# if (is.na(map)) {
				
					# 'hey!'
				
				# # column labels
				# } else if (map %in% c('leftCol', 'middleCol', 'rightCol')) {

					# plot.window(xlim=c(-1, 1), ylim=c(-1, 1))
					# # if (map == 'ul' | substr(map, 1, 1) == 'm') text(0, -0.70, labels=lab, cex=1.4, xpd=NA, pos=1, font=2)
					# # if (substr(map, 2, 2) == 'l') text(0, -0.50, labels=lab, cex=1.4, xpd=NA, pos=1, font=2)
					# text(0, -0.90, labels=lab, cex=1.4, xpd=NA, pos=1, font=2)
					
				# # row labels
				# } else if (map %in% c('topRow', 'middleRow', 'bottomRow')) {

					# plot.window(xlim=c(-1, 1), ylim=c(-1, 1))
					# text(1, 0, labels=lab, cex=1.4, xpd=NA, srt=90, font=2)
					
				# } else {
					
					# rastFileName <- if (is.na(map)) { NA } else
						# if (map == 'ul') { 'defoLocationProb2014to2015' } else
						# if (map == 'um') { 'defoLocationProb2049to2050' } else
						# if (map == 'ur') { 'defoLocationProb2069to2070' } else
						# if (map == 'ml') { 'defoLocationProb2014to2015' } else
						# if (map == 'mm') { 'defoLocationProb2049to2050' } else
						# if (map == 'mr') { 'defoLocationProb2069to2070' } else
						# if (map == 'll') { 'defoLocationProb2014to2015' } else
						# if (map == 'lm') { 'defoLocationProb2049to2050' } else
						# if (map == 'lr') { 'defoLocationProb2069to2070' }
					
					# rastFileName <- paste0('./Deforestation Models/Deforestation Location Probability 2015-2080/', rastFileName, '.tif')
				
					# # get bounding box
					# thisFocus <- paste0('madFocus', inset)

					# focus <- get(thisFocus)
					# ext <- extent(focus)
			
					# # crop geo data
					# pasCrop <- crop(pas, focus)
					# pasCrop <- gUnaryUnion(pasCrop) 
					# humidForestCrop_utm38s <- crop(humidForest_utm38s, focus)
					# humidForestBufferMaskCrop_utm38s <- crop(humidForestBufferMask_utm38s, focus)
					# madagascarCrop_utm38s <- crop(madagascar_utm38s, focus)
					
					# humidForestBufferMaskCrop_utm38s <- calc(humidForestBufferMaskCrop_utm38s, fun=function(x) ifelse(is.na(x), 1, NA))
					
					# x <- raster(rastFileName)
					# xCrop <- crop(x, focus)
					# vals <- cellStats(xCrop, 'sum')
					# insetCols <- if (vals == 0) { cols[1] } else { cols }

					# xMax <- maxValue(x)
					
					# # focal plot
					# plot(humidForestCrop_utm38s, lwd=0.5)
					# plot(humidForestBufferMaskCrop_utm38s, add=TRUE, col='gray80', legend=FALSE)
					# if (!is.null(pasCrop)) plot(pasCrop, border='blue', col=alpha('blue', 0.3), add=TRUE)
					# plot(xCrop, col=cols, breaks=seq(0, xMax, length.out=21), legend=FALSE, add=TRUE)
					# if (!is.null(pasCrop)) plot(pasCrop, border='blue', col=alpha('blue', 0.1), add=TRUE)
					# plot(humidForestCrop_utm38s, lwd=0.5, add=TRUE)
					# plot(madagascarCrop_utm38s, lwd=0.5, border='black', xpd=NA, add=TRUE)
					# plot(focus, border='black', lwd=0.5, add=TRUE)

					# xRange <- ext@xmax - ext@xmin
					# yRange <- ext@ymax - ext@ymin

					# # scale bar
					# if (map %in% c('ul', 'ml', 'll')) {
						
						# scale <- 25000 # meters
						
						# xStart <- 0.05 * xRange + ext@xmin
						# yAt <- 0.03 * yRange + ext@ymin
						
						# lines(c(xStart, xStart + scale), c(yAt, yAt), lwd=4, lend=2, xpd=NA)
						# text(xStart + 0.5 * scale, yAt + 0.05 * yRange, labels='25 km', cex=1.2)
						
					# }

					# # title
					# usr <- par('usr')
					# titleX <- 0.01 * xRange + ext@xmin
					# titleY <- 0.92 * yRange + ext@ymin
					
					# text(titleX, titleY, labels=paste0(letter, ')'), pos=4, cex=1.3, xpd=NA, font=2)

					# # # legend
					# # if (map == 'lr') {
						
						# # # # legend
						# # # labs <- c('0', '0.5', '1')
						# # # # paSwatch <- list(swatchAdjY=c(0, 0.18), col=alpha('blue', 0.3), border='blue', labels='PA')
						# # # paSwatch <- list(swatchAdjY=c(0, 0.18), col=NA, border=alpha('blue', 1), labels='PA')
						# # # legendGrad('bottomright', inset=c(0.2, 0.01), width=0.17, height=0.3, labels=labs, cex=0.83, labAdj=0.5, adjX=c(0.1, 0.7), adjY=c(0.30, 0.75), col=colsLegend, title='Suitability', boxBorder=NA, boxBg=NA, xpd=NA, swatches=list(paSwatch))
						
						# # # legend
						# # labs <- c('0', '0.5', '1')
						# # paSwatch <- list(swatchAdjY=c(0, 0.18), col=NA, border=alpha('blue', 1), labels='PA')
						# # legendGrad('bottomright', inset=c(0.2, 0.01), width=0.17, height=0.3, labels=labs, cex=0.83, labAdj=0.5, adjX=c(0.1, 0.7), adjY=c(0.05, 0.75), col=colsLegend, title='Suitability', boxBorder=NA, boxBg=NA, xpd=NA)
						
					# # }

				# } # if map
					
			# } # next map
			
			# title(main=date(), outer=TRUE, cex.sub=1, line=-2)
			
		# dev.off()
				
		# par(pars)
		
say('#########################################################')
say('### create display maps of forest fragmentation class ###')
say('#########################################################')

	# masks and country borders
	load('./Study Region & Masks/UTM 38S 30-m Resolution/Eastern Humid Forest Polygon.RData')
	load('./Study Region & Masks/UTM 38S 30-m Resolution/Madagascar from GADM 3.6.RData')

	# # # protected areas
	# # pas <- shapefile('./Data/Protected Areas/WDPA_Sept2018_MDG-shapefile-polygons')
	# # pas <- sp::spTransform(pas, CRS(madEaProj))
	# # pas <- crop(pas, madagascar_utm38s)
	# # pas <- gUnaryUnion(pas)
	load('./Data/Protected Areas/WDPA_Sept2018_MDG-shapefile-polygons-onlyTerrestrial.RData')
	
	# graphics parameters
	cols <- c(NA, '#d7191c', '#fdae61', 'yellow', '#a6d96a', 'forestgreen')
	names(cols) <- c('no forest', 'patch', 'transitional', 'perforated', 'edge', 'interior')

	pars <- par()

	# ### Madagascar
	# ##############
	
		# ext <- extent(humidForest_utm38s)
		# ext <- as(ext, 'SpatialPolygons')
		# projection(ext) <- madEaProj
		# ext <- gBuffer(ext, width=100000)
		# pas <- crop(pas, humidForest_utm38s)
		
		# ### setup graticule
		# longs <- seq(43, 51, by=2)
		# lats <- seq(-12, -26, by=-2)
		
		# latLim <- range(lats)
		# longLim <- range(longs)
		
		# grats <- graticule(longs, lats, proj=madEaProj, xlim=longLim, ylim=latLim)
	
		# for (year in c(2014, 2050, 2070)) {
		# # for (year in c(2014)) {
		
			# defos <- if (year <= 2014) { 'anywhere' } else {c('anywhere', 'notPAs') }
		
			# for (defo in defos) {
			# # for (defo in defos[1]) {
		
				# say(year, ' ', defo)
				# outDir <- './Figures & Tables/Forest - Forest Fragmentation Class'
				# dirCreate(outDir)
				
				# if (year <= 2014) {
					# x <- raster(paste0('./Data/Forest - Vieilledent et al 2018/forestFragClass', year, '_utm38s.tif'))
					# pngName <- paste0('Forestation Fragmentation Class for ', year)
				# } else if (defo == 'anywhere') {
					# x <- raster(paste0('./Deforestation Models/Forest Fragmentation Class 2015-2080 Assuming Pop-Sensitive Loss Amount/forestFragClass', year, '_utm38s.tif'))
					# pngName <- paste0('Forestation Fragmentation Class for ', year, ' Assuming Deforestation Anywhere')
				# } else if (defo == 'notPAs') {
					# x <- raster(paste0('./Deforestation Models/Forest Fragmentation Class 2015-2080 Assuming Pop-Sensitive Loss Amount PAs Have 2014 Cover/forestFragClass', year, '_utm38s.tif'))
					# pngName <- paste0('Forestation Fragmentation Class for ', year, ' Assuming PAs have 2014 Cover')
				# }
				
				# png(paste0(outDir, '/', pngName, '.png'), width=4 * 600, height=4 * 800, res=900)

					# par(mar=0.1 * rep(1, 4), oma=0.3 * rep(1, 4))

					# # eastern moist forest
					# plot(ext, border=NA, col=NA)
					# plot(humidForest_utm38s, lwd=0.3, add=TRUE)
					# plot(madagascar_utm38s, lwd=0.3, col='white', add=TRUE)
					# plot(humidForest_utm38s, lwd=0.3, col='white', add=TRUE)
					# plot(x, col=cols, legend=FALSE, add=TRUE)
					# plot(pas, border='blue', lwd=0.3, add=TRUE)
					# plot(grats, lwd=0.3, col='gray', add=TRUE)

					# # insets
					# plot(madFocus1, lwd=0.6, add=TRUE)
					# plot(madFocus2, lwd=0.6, add=TRUE)
					# plot(madFocus3, lwd=0.6, add=TRUE)

					# # graticule labels for longitude
					# longs <- seq(47, 51, by=2)
					# lats <- rep(-12, length(longs))
					
					# longLim <- range(longs)
					# latLim <- range(lats)
					
					# gratLabs <- graticule_labels(longs, lats, xline=min(longLim), yline=max(latLim), proj=madEaProj)
					# longLabs <- gratLabs[gratLabs$islon, ]
					# text(longLabs, label=parse(text=longLabs$lab), col='gray', cex=2.4 * 0.17, srt=90, xpd=NA, adj=c(-0.2, 0.5))
					
					# # graticule labels for latitude
					# lats <- seq(-14, -26, by=-2)
					# longs <- rep(46, length(lats))
					
					# longLim <- range(longs)# + c(-0.02, 0.02)
					# latLim <- range(lats)# + c(-0.02, 0.02)
					
					# gratLabs <- graticule_labels(longs, lats, xline=min(longLim), yline=max(latLim), proj=madEaProj)
					# latLabs <- gratLabs[!gratLabs$islon, ]
					# text(latLabs, label=parse(text=latLabs$lab), col='gray', cex=2.4 * 0.17, xpd=NA, adj=c(1, -0.17))
					
					# # # legend
					# # legend('bottomright', inset=c(0, 0), legend=c(names(cols), 'protected'), fill=c(cols, alpha('blue', 0.35)), border=c(rep('black', length(cols)), 'blue'), cex=0.6, bty='n')
					
				# dev.off()
				
			# } # next defo
			
		# } # next year
		
		# par(pars)

	# ### insets (one plot for all three insets and all three periods)
	# ################################################################

		# humidForestBufferMask_utm38s <- raster('./Study Region & Masks/UTM 38S 30-m Resolution/humidForestBufferMask_utm38s.tif')
	
		# # list of maps
		# maps <- c(NA, 'ul', 'um', 'ur', 'ml', 'mm', 'mr', 'll', 'lm', 'lr', NA, 'leftCol', 'middleCol', 'rightCol', 'topRow', 'middleRow', 'bottomRow')

		# outDir <- './Figures & Tables/Forest - Forest Fragmentation Class'
		# dirCreate(outDir)

		# # future year--left column is always 2014
		# for (futYear in c(2050, 2070)) {
		# # for (futYear in c(2050)) {
			
			# pngName <- paste0('Forest - Forest Fragmentation Class for ', futYear)
			# png(paste0(outDir, '/', pngName, ' Insets.png'), width=3600, height=3600, res=900)
			
				# pars <- par(mfrow=c(4, 4), oma=c(1, 0, 0, 1), mar=c(1, 0, 0, 0))

				# # by EACH RASTER
				# for (map in maps) {
				# # for (map in maps[1:2]) {
					
					# say(futYear, 'subplot ', map)
					
					# lab <- if (is.na(map)) {
						# NA
					# } else if (map == 'topRow') {
						# 'Makira\nNatural Park'
					# } else if (map == 'middleRow') {
						# paste0('Ankeniheny-\nZahamena\nCorridor')
					# } else if (map == 'bottomRow') {
						# paste0('Fandriana-\nVondrozo\nCorridor')
					# } else if (map == 'leftCol') {
						# paste0('2014')
					# } else if (map == 'middleCol') {
						# paste0(futYear, ' forest\nstrict\nprotection')
					# } else if (map == 'rightCol') {
						# paste0(futYear, ' forest\nrelaxed\nprotection')
					# }
					
					# rowCol <- if (is.na(map)) { c(1, 1) } else
						# if (map == 'topRow') { c(2, 1)} else
						# if (map == 'middleRow') { c(3, 1)} else
						# if (map == 'bottomRow') { c(4, 1)} else
						# if (map == 'leftCol') { c(1, 2)} else
						# if (map == 'middleCol') { c(1, 3)} else
						# if (map == 'rightCol') { c(1, 4)} else
						# if (map == 'ul') { c(1, 1) + 1 } else
						# if (map == 'um') { c(1, 2) + 1 } else
						# if (map == 'ur') { c(1, 3) + 1 } else
						# if (map == 'ml') { c(2, 1) + 1 } else
						# if (map == 'mm') { c(2, 2) + 1 } else
						# if (map == 'mr') { c(2, 3) + 1 } else
						# if (map == 'll') { c(3, 1) + 1 } else
						# if (map == 'lm') { c(3, 2) + 1 } else
						# if (map == 'lr') { c(3, 3) + 1 }
					
					# i <- 3
					# letter <- if (is.na(map)) { NA } else
						# if (map == 'ul') { letters[1 + i] } else
						# if (map == 'um') { letters[2 + i] } else
						# if (map == 'ur') { letters[3 + i] } else
						# if (map == 'ml') { letters[4 + i] } else
						# if (map == 'mm') { letters[5 + i] } else
						# if (map == 'mr') { letters[6 + i] } else
						# if (map == 'll') { letters[7 + i] } else
						# if (map == 'lm') { letters[8 + i] } else
						# if (map == 'lr') { letters[9 + i] }
					
					# inset <- if (is.na(map)) { NA } else
						# if (map == 'ul') { 1 } else
						# if (map == 'um') { 1 } else
						# if (map == 'ur') { 1 } else
						# if (map == 'ml') { 2 } else
						# if (map == 'mm') { 2 } else
						# if (map == 'mr') { 2 } else
						# if (map == 'll') { 3 } else
						# if (map == 'lm') { 3 } else
						# if (map == 'lr') { 3 }
					
					# par(mfg=rowCol)
					# plot.new()
					
					# # top left corner
					# if (is.na(map)) {
					
						# 'hey!'
					
					# # column labels
					# } else if (map %in% c('leftCol', 'middleCol', 'rightCol')) {

						# plot.window(xlim=c(-1, 1), ylim=c(-1, 1))
						# text(0, -0.9, labels=lab, cex=1, xpd=NA)
						
					# # row labels
					# } else if (map %in% c('topRow', 'middleRow', 'bottomRow')) {

						# plot.window(xlim=c(-1, 1), ylim=c(-1, 1))
						# text(0.8, 0, labels=lab, cex=1, xpd=NA, srt=90)
						
					# } else {

						# # get raster of interest
						# rastPath <- if (is.na(map)) { NA } else
							# if (map == 'ul') { './Data/Forest - Vieilledent et al 2018' } else
							# if (map == 'um') { './Deforestation Models/Forest Fragmentation Class 2015-2080 Assuming Pop-Sensitive Loss Amount PAs Have 2014 Cover' } else
							# if (map == 'ur') { './Deforestation Models/Forest Fragmentation Class 2015-2080 Assuming Pop-Sensitive Loss Amount' } else
							# if (map == 'ml') { './Data/Forest - Vieilledent et al 2018' } else
							# if (map == 'mm') { './Deforestation Models/Forest Fragmentation Class 2015-2080 Assuming Pop-Sensitive Loss Amount PAs Have 2014 Cover' } else
							# if (map == 'mr') { './Deforestation Models/Forest Fragmentation Class 2015-2080 Assuming Pop-Sensitive Loss Amount' } else
							# if (map == 'll') { './Data/Forest - Vieilledent et al 2018' } else
							# if (map == 'lm') { './Deforestation Models/Forest Fragmentation Class 2015-2080 Assuming Pop-Sensitive Loss Amount PAs Have 2014 Cover' } else
							# if (map == 'lr') { './Deforestation Models/Forest Fragmentation Class 2015-2080 Assuming Pop-Sensitive Loss Amount' }
						
						# rastPathFileName <- if (map %in% c('ul', 'ml', 'll')) {
							# paste0(rastPath, '/forestFragClass2014_utm38s.tif')
						# } else {
							# paste0(rastPath, '/forestFragClass', futYear, '_utm38s.tif')
						# }
					
						# x <- raster(rastPathFileName)
					
						# # get bounding box
						# thisFocus <- paste0('madFocus', inset)

						# focus <- get(thisFocus)
						# ext <- extent(focus)
				
						# # crop geo data
						# pasCrop <- crop(pas, focus)
						# humidForestCrop_utm38s <- crop(humidForest_utm38s, focus)
						# humidForestBufferMaskCrop_utm38s <- crop(humidForestBufferMask_utm38s, focus)
						# madagascarCrop_utm38s <- crop(madagascar_utm38s, focus)
						
						# humidForestBufferMaskCrop_utm38s <- calc(humidForestBufferMaskCrop_utm38s, fun=function(x) ifelse(is.na(x), 1, NA))
						
						# xCrop <- crop(x, focus)
						# vals <- cellStats(xCrop, 'sum')
						# insetCols <- if (vals == 0) { cols[1] } else { cols }

						# xMax <- maxValue(x)
						
						# ### initiate plot
						# plot(humidForestCrop_utm38s, lwd=0.5)

						# ### graticule labels and tick marks: latitude
						# focus_wgs84 <- sp::spTransform(focus, getCRS('wgs84', TRUE))
						# ext <- extent(focus_wgs84)
						# lats <- c(ext@ymin, ext@ymax)
						
						# lats <- round(lats)
						# lats <- seq(lats[1], lats[2], by=0.5)
						# lats <- lats[lats >= ext@ymin & lats <= ext@ymax]
						
						# longs <- rep(ext@xmax, length(lats))
						
						# latLim <- range(lats)
						# longLim <- range(longs)
						
						# if (map %in% c('ur', 'mr', 'lr')) {

							# gratLabs <- graticule_labels(longs, lats, xline=min(longLim), yline=max(latLim), proj=madEaProj)
							# latLabs <- gratLabs[!gratLabs$islon, ]
							# text(latLabs, label=parse(text=latLabs$lab), col='gray', cex=0.6, xpd=NA, adj=c(-0.25, 0.5))
							
						# }

						# # latitude tick marks
						# longLim <- ext@xmax + c(0.03, -0.05)
						# grats <- graticule(longs, lats, proj=madEaProj, xlim=longLim, ylim=c(0, 0))
						# plot(grats, col='gray', add=TRUE)

						# ### graticule labels and tick marks: longitude
						# longs <- c(ext@xmin, ext@xmax)
						
						# longs <- round(longs)
						# longs <- seq(longs[1], longs[2], by=0.5)
						# longs <- longs[longs >= ext@xmin & longs <= ext@xmax]
						
						# lats <- rep(ext@ymax, length(longs))
						
						# latLim <- range(lats)
						# longLim <- range(longs)
						
						# gratLabs <- graticule_labels(longs, lats, xline=min(longLim), yline=max(latLim), proj=madEaProj)
						# longLabs <- gratLabs[gratLabs$islon, ]
						# usr <- par('usr')
						# xx <- c(coordinates(longLabs)[ , 1])
						# yy <- rep(usr[3], length(xx))
						# text(xx, yy, label=parse(text=longLabs$lab), col='gray', cex=0.6, xpd=NA, adj=c(0.5, 1.2))
						
						# # longitude tick marks
						# latLim <- ext@ymin + c(-0.04, 0.05)
						# grats <- graticule(longs, lats, proj=madEaProj, xlim=c(0, 0), ylim=latLim)
						# plot(grats, col='gray', add=TRUE)
						
						# plot(humidForestCrop_utm38s, col='white', lwd=0.5, add=TRUE)
						# plot(humidForestBufferMaskCrop_utm38s, add=TRUE, col='gray80', legend=FALSE)
						# if (!is.null(pasCrop)) plot(pasCrop, border='blue', col=alpha('blue', 0.3), add=TRUE)
						# plot(xCrop, col=cols, legend=FALSE, add=TRUE)
						# if (!is.null(pasCrop)) plot(pasCrop, border='blue', col=alpha('blue', 0.1), add=TRUE)
						# plot(humidForestCrop_utm38s, lwd=0.5, add=TRUE)
						# plot(madagascarCrop_utm38s, lwd=0.5, border='black', xpd=NA, add=TRUE)
						# plot(focus, border='black', lwd=0.5, add=TRUE)

						# ext <- extent(focus)
						# xRange <- ext@xmax - ext@xmin
						# yRange <- ext@ymax - ext@ymin

						# # scale bar
						# if (map %in% c('ul', 'ml', 'll')) {
							
							# scale <- 20000 # meters
							
							# xStart <- 0.14 * xRange + ext@xmin
							# yAt <- 0.05 * yRange + ext@ymin
							
							# lines(c(xStart, xStart + scale), c(yAt, yAt), lwd=3, lend=2, xpd=NA)
							# text(xStart + 0.5 * scale, yAt + 0.08 * yRange, labels=paste(scale / 1000, 'km'), cex=0.8)
							
						# }

						# # # title
						# # usr <- par('usr')
						# # titleX <- -0.05 * xRange + ext@xmin
						# # titleY <- 0.9 * yRange + ext@ymin
						
						# # text(titleX, titleY, labels=paste0(letter, ')'), pos=4, cex=1, xpd=NA, font=2)

					# } # if map
						
				# } # next map
				
				# title(main=date(), outer=TRUE, cex.sub=1, line=-2)
				
			# dev.off()
			
		# } # next year
					
		# par(pars)

	# ### legend
	# ##########
	
		# pngName <- paste0('Forest - Forest Fragmentation Class Legend')
		# png(paste0(outDir, '/', pngName, '.png'), width=800, height=800, res=900)

			# par(mar=rep(0, 4), oma=rep(0, 4), bg=NA)
			# plot(0, 0, col=NA, ann=FALSE, xaxt='n', yaxt='n', main='')
		
			# labs <- c(rev(names(cols)), '', 'protected')
			# legendCols <- c(rev(cols), NA, col=alpha('blue', 0.2))
			# borders <- c(rep(NA, length(cols) - 1), 'gray', NA, 'blue')
		
			# # legend
			# legend(
				# 'left',
				# inset=0.2,
				# legend=labs,
				# fill=legendCols,
				# border=borders,
				# bty='n',
				# cex=0.4
			# )

		# dev.off()
	
	
		
say('DONE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ', date(), level=1, deco='%')
