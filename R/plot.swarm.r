plot.swarm <- function(
	swarm,
	plotEpis=TRUE,
	CRS=swarm$CRS,
	maxRadius=max(swarm$neigh$radius, na.rm=TRUE),
	neighBorder=alpha('blue', 0.2),
	neighLwd=1,
	neighLty='solid',
	neighCol=alpha('cornflowerblue', 0.05),
	epiPch=16,
	epiCol='blue',
	epiBg=NA,
	epiCex=0.5,
	thin=ifelse(sum(swarm$neigh$numNeigh) > 100, 100 / sum(swarm$neigh$numNeigh), 1),
	...
) {
# plot.swarm Plots neighborhoods of a swarm model
#
# ARGUMENTS
# swarm		Object of class `swarm`
# plotEpis	Logical, if TRUE then plots points representing centers of epicenters.
# CRS		Character string represnting Coordinate Reference System to be used for plotting. If not supplied then the CRS in the `swarm` object is used.
# maxRadius	Numeric, maximum size of neighborhood radius to plot.
# thin		Approximate proportion of neighborhoods to plot. If 1, then all neighborhoods are plotted, but thay may be so numerous they cover any underlying plot (e.g., of the study region).
#
# ...		Parameters to pass to plot() and points() functions.  Including "add=TRUE" if adding output to an existing plot.
#
# VALUES
# None.  Produces a map of neighborhoods (on a new or pre-existing graphical device window).  Note that if using
# 
# REQUIRED DEPENDANCIES
# sp, rgeos
#
# OPTIONAL DEPENDANCIES
#
#
# BAUHAUS
# 
#
# EXAMPLE
# FUNCTION()
#
# SOURCE	source('C:/ecology/Drive/r/SDM/SSSDM/plot.swarm.r')
#
# TESTING
#
#
# LICENSE
# This document is copyright ©2017 by Adam B. Smith.  This document is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3, or (at your option) any later version.  This document is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. Copies of the GNU General Public License versions are available at http://www.R-project.org/Licenses/.
#
# AUTHOR	Adam B. Smith | Missouri Botanical Garden, St. Louis, Missouri | adamDOTsmithATmobotDOTorg
# DATE		2017-01
# REVISIONS 

############################
## FUNCTIONS AND PACKAGES ##
############################

####################
## PRE-PROCESSING ##
####################

ellipses <- list(...)

##########
## MAIN ##
##########

### plot neighborhoods
######################

	### establish plot area
	pres <- swarm$pres[ , 1:2]
	pres <- SpatialPoints(pres, CRS(swarm$CRS))
	pres <- sp::spTransform(pres, CRS(swarm$eaCRS))
	ext <- extent(pres)
	maxRadius <- max(swarm$neighStats$radius[!is.na(swarm$neighStats$radius) & !is.infinite(swarm$neighStats$radius)])
	ext <- extent(ext@xmin - maxRadius, ext@xmax + maxRadius, ext@ymin - maxRadius, ext@ymax + maxRadius)
	ext <- as(ext, 'SpatialPolygons')
	projection(ext) <- swarm$eaCRS
	ext <- sp::spTransform(ext, CRS(CRS))
	
	main <- if ('main' %in% names(ellipses)) { ellipses$main } else { '' }
	
	thisAdd <- if ('add' %in% names(ellipses)) { TRUE } else { FALSE }
	
	plot(ext, border=NA, col=NA, main=main, add=thisAdd)

	### convert epicenters to spatial points object, get buffer of set radius, and plot
	justNeigh <- swarm$neighStats[swarm$neighStats$type=='delimited', ]
	justNeigh <- justNeigh[order(justNeigh$radius, decreasing=FALSE), ]
	
	epis <- SpatialPoints(justNeigh[ , c('longitude', 'latitude')], CRS(swarm$CRS))
	epis <- sp::spTransform(epis, CRS(swarm$eaCRS))
	
	# plot neighborhoods
	for (i in seq_along(epis)) {

		if (swarm$neighStats$radius[i] <= maxRadius & runif(1) <= thin) {
	
			neigh <- gBuffer(epis[i], width=swarm$neighStats$radius[i], quadsegs=20)
			neigh <- sp::spTransform(neigh, CRS(CRS))

			plot(neigh, xpd=NA, add=TRUE, border=neighBorder, lwd=neighLwd, lty=neighLty, col=neighCol)
			
		}
		
	}

	### plot epicenters
	if (plotEpis) {
	
		epis <- sp::spTransform(epis, CRS(CRS))
		points(epis[which(swarm$neighStats$radius <= maxRadius)], pch=epiPch, col=epiCol, bg=epiBg, cex=epiCex)
		
	}

#####################
## POST-PROCESSING ##
#####################

}
