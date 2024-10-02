cropToNeighborhood <- function(x, epi, radius, CRS, eaCRS, quadsegs=20) {

	# cropToNeighborhood
	#
	# Crops and masks a Raster* object to a (circular) neighborhood
	
	# x			Raster* object to crop and mask
	# epi		Numeric, 2-values: Coordinates of epicenter
	# radius	Numeric: Neighborhood radius
	# CRS		Character string: Coordinate reference system (i.e., proj4 string)
	# eaCRS		Character string: Equal-area coordinate reference system (i.e., proj4 string)
	# quadsegs	Positive integer: Number of segments used to delineat a quarter section of a buffer

	if (class(epi) == 'numeric') epi <- cbind(epi[1], epi[2])
	
	# convert speicenter to spatial object
	if (class(epi) != 'SpatialPoints' & class(epi) != 'SpatialPointsDataFrame') epi <- SpatialPoints(epi, CRS(CRS))
	epiEa <- sp::spTransform(epi, eaCRS)
	
	# create buffer around epicenter
	buffer <- gBuffer(epiEa, width=radius, quadsegs=quadsegs)
	buffer <- sp::spTransform(buffer, CRS)
	
	# crop and mask
	xNeigh <- crop(x, buffer) # domain stack croppped to neighborhood
	bufferMask <- rasterize(buffer, xNeigh, mask=FALSE)
	xNeigh <- xNeigh * bufferMask
	names(xNeigh) <- names(x)
	
	return(xNeigh)
	
}
	
