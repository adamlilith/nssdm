inspectSwarm <- function(
	swarm,
	prop=c('neighWeight', 'neighWeightScaled', 'neighVsOmni'),
	domain=NULL,
	pause=TRUE,
	...
) {
# inspectSwarm	Explore properties of a swarm object like weighting, coverage, etc.
#
# ARGUMENTS
# swarm		Object of class `swarm`
# prop		'neighWeight' ==> boxplot of effective weight of neighborhoods by radius
#			'neighWeightScaled' ==> barplot of mean effective weight of neighborhoods scaled by mumber of neighborhoods of a given size and their areal coverage (both inside/outside the domain)
# domain	Raster* object used for some aspects of prop.
# pause		Logical, if TRUE then pauses between graphs.
# ...		Arguments to pass to plot() and par()
#
# VALUES
# None (side effect is a plot).
# 
# REQUIRED DEPENDANCIES
# 
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
# SOURCE	source('C:/ecology/Drive/R/SDM/SSSDM/inspectSwarm.r')
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

	par(ask=pause)

##########
## MAIN ##
##########

	### neighborhood weights
	########################
	
	if ('neighWeight' %in% prop) {
	
		par(mfrow=c(1, 1))
		
		y <- swarm$neighStats$neighWeight[!is.na(swarm$neighStats$radius)]
		x <- swarm$neighStats$radius[!is.na(swarm$neighStats$radius)]
		x <- as.factor(x)
		
		boxplot(y ~ x, col=alpha('cornflowerblue', 0.5), xlab='Radius', ylab='Approximate weight', main='Effective weight per neighborhood')
	
	}

	if ('neighWeightScaled' %in% prop) {
	
		# calculate mean weights and number of neighborhoods of each size plus neihgborhood area
		neighWeight <- aggregate(swarm$neighStats[ , c('radius', 'neighWeight')], by=list(swarm$neighStats$radius), FUN=sum)$neighWeight
		radii <- unique(na.omit(swarm$neighStats$radius))
		aggArea <- pi * radii^2
		
		scaledWeight <- neighWeight * aggArea
	
		# gaussian weighting
		if (any(swarm$neighStats$weighting == 'gauss')) {
			
			gaussWeight <- dnorm(radii / 2, mean=0, sd=unique(na.omit(swarm$neighStats$sd)))
			scaledWeight <- scaledWeight * gaussWeight
			
		}
	
		# include omnibus models
		if ('omnibus' %in% swarm$neighStats$type) {
		
			omniWeight <- sum(swarm$neighStats$neighWeight[swarm$neighStats$type=='omnibus'])
			
			areaRast <- area(domain[[1]], na.rm=T)
			areaDomain <- cellStats(areaRast, 'sum')
			
			omniWeight <- omniWeight * areaDomain
			
			radii <- c(radii, 'Inf')
			scaledWeight <- c(scaledWeight, omniWeight)
		
		}
		
		par(mfrow=c(1, 2))

		# effective weight per radius
		barplot(scaledWeight, names.arg=radii, xlab='Radius', ylab='Approxmiate total weight', main='Total Weight')
		
		# cumulative weight per radius
		cumulWeight <- cumsum(scaledWeight)
		cumulWeight <- cumulWeight / max(cumulWeight)
		
		barplot(cumulWeight, names.arg=radii, xlab='Radius', ylab='Approximate cumulative proportion of weight', main='Cumulative weight')
		
		title(main='Effective weight scaled by number of neighborhoods and area (and possibly mean Gaussian effect)', outer=TRUE, line=-1, xpd=NA)
		
	}
	
#####################
## POST-PROCESSING ##
#####################

par(ask=FALSE)

}
