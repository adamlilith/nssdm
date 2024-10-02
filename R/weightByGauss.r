weighByGauss <- function(x, epi, sd, distFx, ...) {

	## Returns raster with Gaussian weighting.  For internal use.
	
	# x			raster object
	# epi		SpatialPoints* object representing neighborhood epicenter
	# sd		Standard deviation of Gaussian curve
	# distFx	distance function
	# ...		Extra stuff (not used)

	x <- x * 0
	
	cells <- 1:ncell(x)
	
	xy <- xyFromCell(x, cells, spatial=TRUE)
	
	cellDist <- distFx(epi, xy)
	gaussWeights <- dnorm(cellDist, mean=0, sd=sd)
	
	x[] <- gaussWeights
	
	x
	
}

