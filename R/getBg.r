getBg <- function(
	swarm,
	domain,
	bias=NULL,
	numBg=10000,
	verbose=TRUE,
	...
) {
# getBg		Obtains background points and extracts environmental variables for them.
#
# ARGUMENTS
# swarm		Object of class `swarm`
# domain	RasterStack or RasterBrick of predictors.  Assumed to be in WGS84 (unprojected) datum.
# bias		Object of class `Raster` with cell values equal to relative probability of being drawn.  If NULL then all cells will have equal probability of being drawn (after correcting for cell size).
# numBg		Integer, number of background sites to choose per neighborhood.
# verbose	Logical, if TRUE then displays progress.
#
# VALUES
# Object of class `swarm` with background sites for each neighborhood and associated environmental variables for each site.
# 
# REQUIRED DEPENDANCIES
# sp, raster
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
# SOURCE	source('C:/ecology/Drive/R/SDM/SSSDM/getBg.r')
#
# TESTING
#
#
# LICENSE
# This document is copyright ©2014 by Adam B. Smith.  This document is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3, or (at your option) any later version.  This document is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. Copies of the GNU General Public License versions are available at http://www.R-project.org/Licenses/.
#
# AUTHOR	Adam B. Smith | Missouri Botanical Garden, St. Louis, Missouri | adamDOTsmithATmobotDOTorg
# DATE		2016-12
# REVISIONS 

############################
## FUNCTIONS AND PACKAGES ##
############################

####################
## PRE-PROCESSING ##
####################

# get bias raster if not supplied
if (is.null(bias)) bias <- area(domain[[1]], na.rm=TRUE)

# remember number of BG sites per neighborhood
swarm$minContrast <- numBg

swarm$neighStats$numContrast <- numBg

##########
## MAIN ##
##########

if (verbose) {
	if (verbose) say('=== Getting Background Sites ===')
	pb <- txtProgressBar(min=0, max=nrow(swarm$neighStats), width=30, style=3)
}

for (countNeigh in nrow(swarm$neighStats):1) {

	# if (verbose) say('Getting background sites for neighborhood of radius ', swarm$neighStats$radius[countNeigh], ' (', round(100 * (nrow(swarm$neighStats) - countNeigh) / nrow(swarm$neighStats), 1), '% done) at ', date(), '.')

	### delimited neighborhood
	##########################

	if (swarm$neighStats$type[countNeigh]=='delimited') {
	
		# make buffer around neighborhood epicenter
		thisBias <- cropToNeighborhood(x=bias, epi=swarm$neighStats[countNeigh, c('longitude', 'latitude')], radius=swarm$neighStats$radius[countNeigh], CRS=swarm$CRS, eaCRS=swarm$eaCRS, quadsegs=20)
		
		# throw error if neighborhood has just one cell
		if (cellStats(thisBias * 0 + 1, 'sum')==1) {
			stop('Neighborhood #', countNeigh, ' (radius = ', swarm$neighStats$radius[countNeigh], ') is smaller than a single cell. Models will be not trainable. Increase neighborhood size.')
		}
		
	### omnibus
	###########
	
	} else {
	
		# get background points from entire domain
		thisBias <- bias
	
	}

	# extract background sites
	bg <- sampleRast(mask=thisBias, n=numBg, replace=TRUE, prob=TRUE)
		
	# extract environment at BG sites
	env <- as.data.frame(extract(domain, bg))

	# remember: get index of these BG points
	contrastRows <- if (is.null(swarm$contrast)) {
		1:nrow(bg)
	} else {
		nrow(swarm$contrast) + 1:nrow(bg)
	}
	
	# remember: calculate weights
	globalContrastWeight <- rep(1:nrow(bg))
	
	if (swarm$neighStats$weighting[countNeigh]=='gauss') {
	
		epi <- cbind(swarm$neighStats$longitude[countNeigh], swarm$neighStats$latitude[countNeigh])
		epi <- SpatialPoints(epi, CRS(swarm$CRS))
		bgSp <- SpatialPoints(cbind(bg[ , 1:2]), CRS(swarm$CRS))
		distContrast <- swarm$neigh$distFx(p1=epi, p2=bgSp)
		globalContrastWeight <- globalContrastWeight * dnorm(distContrast, mean=0, sd=swarm$neighStats$sd[countNeigh])
	
	}	
		
	# remember: assign in/out of bag samples
	swarm$neighPresContrast[[countNeigh]]$contrastIndex <- contrastRows
	swarm$neighPresContrast[[countNeigh]]$globalContrastWeight <- globalContrastWeight
	
	inBag <- sort(sample(contrastRows, round(swarm$bag * length(contrastRows))))
	ooBag <- contrastRows[!(contrastRows %in% inBag)]
	
	swarm$neighPresContrast[[countNeigh]]$contrastIndexInBag <- inBag
	swarm$neighPresContrast[[countNeigh]]$contrastIndexOutOfBag <- ooBag

	bg <- as.data.frame(bg)
	names(bg) <- names(swarm$pres)[1:2]
	bg$coverage <- NA
	bg <- cbind(bg, env)
	
	if (is.null(swarm$contrast)) {
		swarm$contrast <- bg
	} else {
		swarm$contrast <- rbind(swarm$contrast, bg)
	}
	
	if (verbose) setTxtProgressBar(pb, nrow(swarm$neighStats) - countNeigh)
	
} # next neighborhood

#####################
## POST-PROCESSING ##
#####################

return(swarm)

}
