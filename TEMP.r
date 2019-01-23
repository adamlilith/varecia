# source('C:/Ecology/Drive/Research/Varecia (Toni Lyn Morelli)/Versions 06/Code/TEMP.r')

	png(paste0(outDir, '/Occurrence Map for Varecia - Both Species.png'), width=600, height=1200, res=600)
		
		par(mar=0.5 * c(1, 1, 1, 1), oma=rep(0, 4))
		
		plot(humidForest_utm38s, ann=FALSE)
		
		for (i in seq_along(grays)) grays[i] <- alpha(grays[i], 0.6)
		plot(hs, add=TRUE, col=grays, legend=FALSE)
		plot(madagascar_utm38s, add=TRUE, lwd=0.7)
		plot(humidForest_utm38s, add=TRUE, border=alpha('chartreuse', 0.4), col=alpha('chartreuse', 0.1), lwd=0.6)
		plot(humidForestBuffer_utm38s, add=TRUE, border='darkred', col=NA, lwd=0.6)
		plot(pas, add=TRUE, border='blue', col=alpha('blue', 0.3), lwd=0.6)
		points(varecia, pch=1, bg=bg, cex=0.2)
		
		legend('bottomright', bty='n', legend=c('V. variegata', 'V. rubra', 'Humid forest', 'Study region', 'Protected'), pch=c(21, 24, NA, NA, NA), col=c('black', 'black', NA, NA, NA), pt.bg=c('white', 'red', NA, NA, NA), border=c(NA, NA, 'chartreuse', 'darkred', 'blue'), fill=c(NA, NA, 'darkseagreen1', NA), cex=0.3)

		title(sub=date(), cex.sub=0.1, line=-0, xpd=NA)
		
	dev.off()
