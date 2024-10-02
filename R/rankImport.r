rankImport <- function(
	x,
	nth=1,
	...
) {
# rankImport For a stack of rasters representing one measure of variable importance (e.g., 'importR', 'importAuc', etc.), returns a single raster with the n-th most important variable indicated.
#
# ARGUMENTS
# x				RasterStack or RasterBrick, one per predictor in a swarm model.
# nth			Rank of variable desired.  If this is 1, then the most important variable will be chosen.  If 2, then the second most important, and so on. In case of ties the first variable in x is used (if ties seem important to the result, it may be helpful to recalculate importance rank using x with variables sorted in a different order).
# ...			Arguments to pass to "calc()" (e.g., `filename`, `format`, `overwrite`, etc.)
#
# VALUES
# Raster object with integer values from 1 to <number of layers in x>, each corresponding to a predictor (in the same order as listed in x).  For example, an output raster might have values of 1 in one set of cells, 2 in another,  etc.  This indicates that the n-th most important variable in the cells with 1s is the variable represented by the first layer in x, the n-th most important variable in cells with 2s is the variable represented by the second layer in x, etc.
# 
# REQUIRED DEPENDANCIES
# raster
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
# SOURCE	source('C:/ecology/Drive/R/SDM/SSSDM/rankImport.r')
#
# TESTING
#
#
# LICENSE
# This document is copyright ©2016 by Adam B. Smith.  This document is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3, or (at your option) any later version.  This document is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. Copies of the GNU General Public License versions are available at http://www.R-project.org/Licenses/.
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

ranked <- calc(x, fun=function(X, na.last='keep', ...) rank(X, ties.method='min')[nth], na.rm=TRUE)
names(ranked) <- paste0('indexOf', nth, 'MostImportVar')

#####################
## POST-PROCESSING ##
#####################

return(ranked)

}

