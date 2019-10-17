# source('C:/Ecology/Drive/Research/Varecia (Toni Lyn Morelli)/Versions 06/Code/TEMP.r')

		# all of Madagascar
		png(paste0(outDir, '/', pngName, '.png'), width=1800, height=1800, res=450)

			par(mfrow=c(1, 2), mar=rep(0, 4), oma=c(0, 0, 1, 0), bg='black')
		
			# for (protection in c('anywhere', 'notPAs')) {
			for (protection in c('anywhere', 'anywhere')) {
			
				say(year, ' ', protection)

				x1 <- raster(paste0('./Deforestation Models/Forest 2015-2080 Assuming Pop-Sensitive Loss Amount/forest', year, '.tif'))
				
				if (protection == 'anywhere') {
					
					label <- 'Relaxed\nprotection'
					
				} else if (protection == 'notPAs') {

					label <- 'Strict\nprotection'

					# add forest in PAs
					x1pas <- stack(x1, forest2014inWpdaPAs)
					x1pas <- max(x1pas, na.rm=TRUE)
						
					fx <- function(x) ifelse(x == 1, 1, NA)
					x1 <- clusterR(x1pas, calc, args=list(fun=fx))
					
				}

				# eastern moist forest
				plot(madagascar_utm38s, border='white', ann=FALSE)
				plot(hs, col=hsCols, legend=FALSE, ann=FALSE, add=TRUE)
				plot(x1, col=x1col, legend=FALSE, add=TRUE)
				plot(notHumidForest_utm38s, col=alpha('black', 0.4), border=NA, add=TRUE, ann=FALSE)
				if (protection == 'notPAs') plot(pas, border=paBorder, lwd=0.4, add=TRUE, ann=FALSE)

				text(450000, 8473849, labels=label, cex=0.8, xpd=NA, pos=4, col='white')
				
			} # next forest protection
			
			title(main=year, cex.main=1, outer=TRUE, line=-11, col.main='white')

		dev.off()
