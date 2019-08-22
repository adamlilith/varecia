# source('C:/ecology/Drive/Research/Varecia (Toni Lyn Morelli)/Versions 06/Code/TEMP.r')

# say('#############################################################')
# say('### create displays of ecological niche model predictions ###')
# say('#############################################################')

	# generalization
	gcm <- 'EnsembleMean'
		
	outDir <- './Figures & Tables/Ecological Niche Models - Predictions/'
	dirCreate(outDir)

	cexMult <- 0.45 # multiplier for cex

	# RCPs for each column
	rcpCol2 <- '2pt6'
	rcpCol3 <- '4pt5'
	rcpCol4 <- '6pt0'
	rcpCol5 <- '8pt5'
	
	rcpNice <- function(rcp) gsub(rcp, pattern='pt', replacement='.')
	rcpCol2Nice <- rcpNice(rcpCol2)
	rcpCol3Nice <- rcpNice(rcpCol3)
	rcpCol4Nice <- rcpNice(rcpCol4)
	rcpCol5Nice <- rcpNice(rcpCol5)
		
	# ancillary geo data
	humidForestBufferMask_utm38s <- raster('./Study Region & Masks/UTM 38S 30-m Resolution/humidForestBufferMask_utm38s.tif')

	load('./Study Region & Masks/UTM 38S 30-m Resolution/Eastern Humid Forest Polygon.RData')
	load('./Study Region & Masks/UTM 38S 30-m Resolution/Madagascar from GADM 3.6.RData')
	
	# # pas <- shapefile('./Data/Protected Areas/WDPA_Sept2018_MDG-shapefile-polygons')
	# # pas <- sp::spTransform(pas, CRS(madEaProj))
	# # pas <- crop(pas, madagascar_utm38s)
	
	load('./Data/Protected Areas/WDPA_Sept2018_MDG-shapefile-polygons-onlyTerrestrial.RData')
	
	# list of maps
	maps <- c('r1c1', 'r2c1', 'r3c1', 'r1c2', 'r2c2', 'r3c2', 'r1c3', 'r2c3', 'r3c3', 'r1c4', 'r2c4', 'r3c4', 'r1c5', 'r2c5', 'r3c5')
	# maps <- c('r1c1', 'r2c1')
	# maps <- c('r3c5')
	
	# color ramp for focal rasters
	cols <- colorRampPalette(c('lightsalmon', 'red', 'darkred'))
	colBreaks <- 5 # number of color categories plus 1
	cols <- cols(colBreaks - 2)
	cols <- c(NA, cols)

	colsLegend <- cols

	# hillshade
	hs <- raster('./Data/Topography - GMTED2010/hillshadeGmted2010_utm38s.tif')
	
	# hillshade colors
	hsCols <- paste0('gray', 0:100)
	hsCols <- alpha(hsCols, 0.7)

	for (futYear in c(2050, 2070)) {
	# for (futYear in c(2050)) {
	
		say(futYear)

		# prediction rasters
		# rasters (u = upper, m = middle, l = lower or left, r = right)
		rastDir <- paste0('./Ecological Niche Models/Prediction Rasters')
		
		r1c1 <- raster(paste0(rastDir, '/glmEnm_vareciaGenus_climateCurrent_forest2014_utm38s.tif'))
		r2c1 <- raster(paste0(rastDir, '/glmEnm_vareciaGenus_climateCurrent_forest', futYear, '_defoOutsidePAs_utm38s.tif'))
		r3c1 <- raster(paste0(rastDir, '/glmEnm_vareciaGenus_climateCurrent_forest', futYear, '_defoAnywhere_utm38s.tif'))

		r1c2 <- raster(paste0(rastDir, '/glmEnm_vareciaGenus_climate', futYear, gcm, 'Rcp', rcpCol2, '_forest2014_utm38s.tif'))
		r2c2 <- raster(paste0(rastDir, '/glmEnm_vareciaGenus_climate', futYear, gcm, 'Rcp', rcpCol2, '_forest', futYear, '_defoOutsidePAs_utm38s.tif'))
		r3c2 <- raster(paste0(rastDir, '/glmEnm_vareciaGenus_climate', futYear, gcm, 'Rcp', rcpCol2, '_forest', futYear, '_defoAnywhere_utm38s.tif'))

		r1c3 <- raster(paste0(rastDir, '/glmEnm_vareciaGenus_climate', futYear, gcm, 'Rcp', rcpCol3, '_forest2014_utm38s.tif'))
		r2c3 <- raster(paste0(rastDir, '/glmEnm_vareciaGenus_climate', futYear, gcm, 'Rcp', rcpCol3, '_forest', futYear, '_defoOutsidePAs_utm38s.tif'))
		r3c3 <- raster(paste0(rastDir, '/glmEnm_vareciaGenus_climate', futYear, gcm, 'Rcp', rcpCol3, '_forest', futYear, '_defoAnywhere_utm38s.tif'))

		r1c4 <- raster(paste0(rastDir, '/glmEnm_vareciaGenus_climate', futYear, gcm, 'Rcp', rcpCol4, '_forest2014_utm38s.tif'))
		r2c4 <- raster(paste0(rastDir, '/glmEnm_vareciaGenus_climate', futYear, gcm, 'Rcp', rcpCol4, '_forest', futYear, '_defoOutsidePAs_utm38s.tif'))
		r3c4 <- raster(paste0(rastDir, '/glmEnm_vareciaGenus_climate', futYear, gcm, 'Rcp', rcpCol4, '_forest', futYear, '_defoAnywhere_utm38s.tif'))

		r1c5 <- raster(paste0(rastDir, '/glmEnm_vareciaGenus_climate', futYear, gcm, 'Rcp', rcpCol5, '_forest2014_utm38s.tif'))
		r2c5 <- raster(paste0(rastDir, '/glmEnm_vareciaGenus_climate', futYear, gcm, 'Rcp', rcpCol5, '_forest', futYear, '_defoOutsidePAs_utm38s.tif'))
		r3c5 <- raster(paste0(rastDir, '/glmEnm_vareciaGenus_climate', futYear, gcm, 'Rcp', rcpCol5, '_forest', futYear, '_defoAnywhere_utm38s.tif'))

		png(paste0(outDir, '/ENM Predictions Using Ensemble Mean across GCMs for ', futYear, ' x Deforestation x RCPs.png'), width=4 * 400, height=3 * 900, res=900)
		
			pars <- par(mfrow=c(3, 5), oma=c(0, 1, 2, 0), mar=c(0, 0, 0, 0))

			for (map in maps) {
				
				say('subplot ', map)
				
				rowCol <- if (map == 'r1c1') { c(1, 1)} else
					if (map == 'r1c2') { c(1, 2)} else
					if (map == 'r1c3') { c(1, 3)} else
					if (map == 'r1c4') { c(1, 4)} else
					if (map == 'r1c5') { c(1, 5)} else
				
					if (map == 'r2c1') { c(2, 1)} else
					if (map == 'r2c2') { c(2, 2)} else
					if (map == 'r2c3') { c(2, 3)} else
					if (map == 'r2c4') { c(2, 4)} else
					if (map == 'r2c5') { c(2, 5)} else
				
					if (map == 'r3c1') { c(3, 1)} else
					if (map == 'r3c2') { c(3, 2)} else
					if (map == 'r3c3') { c(3, 3)} else
					if (map == 'r3c4') { c(3, 4)} else
					if (map == 'r3c5') { c(3, 5)}
				
				letter <- if (is.na(map)) { NA } else
					if (map == 'r1c1') { 'a' } else
					if (map == 'r1c2') { 'b' } else
					if (map == 'r1c3') { 'c' } else
					if (map == 'r1c4') { 'd' } else
					if (map == 'r1c5') { 'e' } else

					if (map == 'r2c1') { 'f' } else
					if (map == 'r2c2') { 'g' } else
					if (map == 'r2c3') { 'h' } else
					if (map == 'r2c4') { 'i' } else
					if (map == 'r2c5') { 'j' } else
				
					if (map == 'r3c1') { 'k' } else
					if (map == 'r3c2') { 'l' } else
					if (map == 'r3c3') { 'm' } else
					if (map == 'r3c4') { 'n' } else
					if (map == 'r3c5') { 'o' }
				
				par(mfg=rowCol)
				plot.new()
				
				x <- get(map)
				
				# main plot
				plot(humidForest_utm38s, col=NA, border='black', lty='solid', lwd=0.2)
				plot(pas, border=NA, col='lightskyblue', add=TRUE)
				# # plot(pas, border=alpha('blue', 0.4), col=alpha('blue', 0.175), lwd=0.4, add=TRUE)
				plot(x, col=cols, breaks=seq(0, 100, length.out=colBreaks), legend=FALSE, add=TRUE)
				plot(pas, border='lightskyblue', lwd=0.3, add=TRUE)
				plot(madagascar_utm38s, lwd=0.4, add=TRUE)

				# title
				usr <- par('usr')
				titleX <- usr[1] - 0 * (usr[2] - usr[1])
				titleY <- usr[4] - 0.1 * (usr[4] - usr[3])
				text(titleX, titleY, labels=paste0(letter, ')'), pos=4, cex=cexMult * 1.1, xpd=NA)

				# legend
				if (map == 'r3c5') {

					# legend
					labs <- c(0.25, 0.5, 0.75)
					paSwatch <- list(swatchAdjY=c(0, 0.13), col='lightskyblue', border='black', labels='PA', lwd=0.2)

					legendBreaks('bottomright', inset=c(0.28, 0.005), width=0.15, height=0.25, labels=labs, labAdjX=-0.3, labAdjY=labs, cex=cexMult * 0.73, adjX=c(0.1, 0.7), adjY=c(0.22, 0.75), col=cols, colBorder=NA, title='Suitability', boxBorder=NA, boxBg=NA, xpd=NA, swatches=list(paSwatch), lwd=0.2)
					
				}

				### insets
				
				insetLwd <- 0.5
				plot(madFocus1, lwd=insetLwd, lend=1, add=TRUE)
				plot(madFocus2, lwd=insetLwd, lend=1, add=TRUE)
				plot(madFocus3, lwd=insetLwd, lend=1, add=TRUE)

				# column label (RCPs)
				if (map %in% c('r1c1', 'r1c2', 'r1c3', 'r1c4', 'r1c5')) {

					labCol <- if (is.na(map)) {
						NA
					} else if (map == 'r1c1') {
						'Current climate'
					} else if (map == 'r1c2') {
						paste0(futYear, ' climate\nRCP ', rcpCol2Nice)
					} else if (map == 'r1c3') {
						paste0(futYear, ' climate\nRCP ', rcpCol3Nice)
					} else if (map == 'r1c4') {
						paste0(futYear, ' climate\nRCP ', rcpCol4Nice)
					} else if (map == 'r1c5') {
						paste0(futYear, ' climate\nRCP ', rcpCol5Nice)
					}
					
					x <- 0.5 * (usr[1] + usr[2])
					y <- usr[4] + 0.12 * (usr[4] - usr[3])
					text(x, y, labels=labCol, cex=cexMult * 1, xpd=NA)
				
				}
				
				# row labels (defo)
				if (map %in% c('r1c1', 'r2c1', 'r3c1')) {

				usr <- par('usr')
					labRow <- if (map == 'r1c1') {
						'2014 forest'
					} else if (map == 'r2c1') {
						paste0(futYear, ' forest\nstrict protection')
					} else if (map == 'r3c1') {
						paste0(futYear, ' forest\nrelaxed protection')
					}
				
					x <- usr[1] - 0.22 * (usr[2] - usr[1])
					y <- 0.5 * (usr[3] + usr[4])
					text(x, y, labels=labRow, cex=cexMult * 1.1, xpd=NA, srt=90)
					
				}
					
			} # next map
			
		dev.off()

	} # next future year

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
	
