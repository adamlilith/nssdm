makeSwarm <- function(
	domain,
	eaCRS,
	pres,
	contrast=NULL,
	bag=0.7,
	minPres=30,
	minContrast=NULL,
	minPresWeight=NULL,
	minContrastWeight=NULL,
	rawPresWeights=NULL,
	rawContrastWeights=NULL,
	neigh=list(
		weighting='gauss',
		radius=1000 * 2^(1:10),
		sd=0.5 * 1000 * 2^(1:10),
		numDesired=2^(10:1),
		rawNeighWeights=2^(10:1),
		epiType='random',
		attempts=100,
		omnibus=1,
		omnibusWeightScale=2^-1,
		distFx=distCosine
	),
	verbose=TRUE
) {
# makeNeigh Produces a template object of class "swarm" with a "neighborhood" data frame that contains coordinates for epicenters and radii plus neighborhood weights.
#
# ARGUMENTS (aside from those already in sdsdm())
#
# domain		Raster stack of predictors.  Assumed to be in WGS84 (unprojected) datum.
# eaCRS			proj4 string for an equal-area projection relevant to the domain (e.g., Albers Equal-Area North American)
# modelDir		directory in which to save models
# pres			Data frame or matrix with presence sites used to train model, 1st column is longitude and 2nd is latitude.  If there are more than two columns, the remainder are assumed to contain environmental data pertaining to each site.
# contrast		NULL OR data frame or matrix with absence/background sites used to train model, in which case the 1st column is longitude and 2nd is latitude.  If there are more than two columns, the remainder are assumed to contain environmental data pertaining to each site.
# bag			Fraction of presences (and absences if contrast is not NULL) within a neighborhood to use to train a submodel.
# minPres		Minimum number of presences (*before* removing 1 - bag portion) needed within a neighborhood.
# minContrast		Minimum number of absences (*before* removing 1 - bag portion if absences are supplied by user) needed within a neighborhood to train a model.  Ignored if contrast = NULL.  Note that this has the potential to limit the lower size of a viable neighborhood because if no neighborhood can contain minContrast absence sites, then the size of the neighborhood is too small.
# minPresWeight Minimum total weight of presences in a neighborhood to train a model (*before* removing 1 - bag portion). If NULL then this requirement for a neighborhood is ignored.
# minContrastWeight Minimum total weight of absences (if supplied) in a neighborhood to train a model (*before* removing 1 - bag portion). If NULL then this is equal to minContrast.  Ignored if NULL.
# rawPresWeights	Numeric vector same length as nrow(pres) with weights of presences.  If NULL each presence is assumed to have equal weight (before application of neighborhood weighting, if any).
# rawContrastWeights	Numeric vector same length as nrow(contrast) with weights of absences.  If NULL each absence is assumed to have equal weight (before application of any neighborhood weighting, if any).
# neigh			List item used to define neighborhood around each epicenter.  Contains different elements, depending on type of neighborhood.  Examples:
#				  Each presences/contrast point receives a weight of 1 and neihgborhood radii are 50 km: neigh=list(weighting='none', radius=50000)
#				  Draw radii with probability proportional to a user-supplied set of weights: neigh=list(weighting='none', radius=c(50000, 100000, 200000), drawBy=c(3, 2, 1))
#				  Weighting of *all* presences and contrast sites up to 60000 km according to Gaussian kernel with sd=20 km: neigh=list(weighting='gauss', sd=20, radius=60000)
#				  Draw randomly from a set of distances with a probability drawn proportionate to a user-defined set of weights: neigh=list(weighting='gauss', radius=c(10000, 20000, 50000), sd=c(5000, 10000, 25000), drawBy=c(2, 1, 1))
#
# NOTE: If there is more than radius value supplied (length(radius) > 1), the script will ONLY allow randomly located epicenters and will attempt to place the epicenter until either there are enough points in the neighborhood to create a model or until "attempts" epicenters have been located (and failed).  If the latter happens, a new radius is chosen (though it may possibly be the same size as the last).
#
# epicenters	Either an integer indicating number of epicenters to randomly locate or 2-column data frame or matrix with coordinates for epicenters.
# attempts		Used only if epicenters are randomly located or if length(neigh$radius) > 1.  Number of times to attempt either to locate an epicenter or select a neighborhood radius that has sufficient training presences and contrast points.
# verbose		TRUE ==> Display detailed progress.  FALSE ==> Display minimal progress.
#
# epiType
# 'random' ==> randomly place epicenters number of epicenters, each of which has minPres number of presences (after out-of-bag fraction is removed).  Note that this does not ensure all presences will be covered by at least one neighborhood.
# 'inclusion' ==> place neighborhoods so that each presence is contained within at least minCover neighborhoods.
#
# stopRule
# 'cover' ==> stop when each presence (and absence) is contained in at least minCover neighborhoods
# 'number' ==> stop when there are numEpi neighborhoods that have been successfully placed
#
# minCover
# Used only if stopRule=='cover'.  Integer, number of times each presence (and absence, if supplied) should be included in a neighborhood.  Note that multiple attempts will be made to include a site, but in some cases, depending on the spatial distribution of sites and radii chosen, it will not be possible to include a site.
#
# numEpi
# Used only if stopRule=='number'.  Integer, number of epicenters to place.  Note that this will not ensure (or try to ensure) all presences (and absences) are in at least one neighborhood.
#
# BAUHAUS
# For random neighborhoods, locate epicenter, choose radius, see if it has enough presences, if not increase radius, choose another if not possible
# For coverage neighborhoods, locate presence with < minCover neighborhoods, locate epicenter within a given radius, see if has enough presences, increase radius if not (to max, then give warning)
#
# VALUES
# an object of class 'swarm'

############################
## FUNCTIONS AND PACKAGES ##
############################

###############################
### WARNINGS/ERROR CATCHING ###
########################### ###

# flag to indicate if absences are supplied by user or will be background sites (chosen by training script)
contrastSupplied <- !is.null(contrast)

# catch input errors
if (contrastSupplied) {
	
	if (minContrast > nrow(contrast)) stop('Minimum number of required contrast sites (minContrast) is > number of supplied contrast sites.')
	if (minContrast > bag * nrow(contrast)) stop('Minimum number of required contrast sites (minContrast) is > number of supplied contrast sites (nrow(contrast) * bag).')

}

if (is.na(projection(domain))) stop('Domain rasters are lacking a coordinate reference system (CRS).')

####################
## PRE-PROCESSING ##
####################

# convert presenses (and absences) to data frames
if (class(pres)=='matrix') pres <- as.data.frame(pres)
if (contrastSupplied) if (class(contrast)=='matrix') contrast <- as.data.frame(contrast)

# if not supplied, generate a priori weights for presences and absences
if (is.null(rawPresWeights)) rawPresWeights <- rep(1, nrow(pres))
if (contrastSupplied & is.null(rawContrastWeights)) rawContrastWeights <- rep(1, nrow(contrast))

# if not supplied, generate minimum weights
if (is.null(minPresWeight)) minPresWeight <- 0
if (contrastSupplied & is.null(minContrastWeight)) minContrastWeight <- 0

# add auxillary variables to presences (and absences)
pres$coverage <- 0 # number of times each presence is included in a neighborhood

if (contrastSupplied) {
	contrast$coverage <- 0 # number of times each absence is included in a neighborhood
}

# get environmental data at presences (and absences)
if (ncol(pres)==3) pres <- cbind(pres, as.data.frame(extract(domain, pres[ , 1:2])))
if (contrastSupplied) if (ncol(contrast)==3) contrast <- cbind(contrast, as.data.frame(extract(domain, contrast[ , 1:2])))

# get equal-area versions of presences (and absences)
presSp <- SpatialPoints(pres[ , 1:2], CRS(projection(domain)))
presEa <- sp::spTransform(presSp, CRS(eaCRS))
if (contrastSupplied) {
	contrastSp <- SpatialPoints(pres[ , 1:2], CRS(projection(domain)))
	contrastEa::spTransform(contrastSp, CRS(eaCRS))
}

# initiate swarm object... use "NA" for placeholder of empty items
swarm <- list()
class(swarm) <- c('swarm', 'list')
swarm$start <- NA
swarm$end <- NA
swarm$domainNames <- names(domain)
swarm$CRS <- projection(domain)
swarm$eaCRS <- eaCRS
swarm$eaExtent <- extent(domain)
swarm$modelDir <- NA
swarm$pres <- pres
swarm$contrast <- contrast
swarm$rawPresWeights <- rawPresWeights
swarm$rawContrastWeights <- if (!is.null(rawContrastWeights)) { rawContrastWeights } else { NA }
swarm$bag <- bag
swarm$minPres <- minPres
swarm$minContrast <- minContrast
swarm$minPresWeight <- minPresWeight
swarm$minContrastWeight <- minContrastWeight
swarm$bias <- NA
swarm$model <- NA
swarm$neigh <- neigh
swarm$neigh$numNeigh <- rep(0, length(neigh$radius))
swarm$neighStats <- data.frame(type=character(), weighting=character(), radius=numeric(), sd=numeric(), neighWeight=numeric(), numPres=integer(), numContrast=integer(), longitude=numeric(), latitude=numeric())
swarm$neighPresContrast <- list()

# if using random epicenter location, create buffer around presences (and absences) in which to locate epicenters... enables selection of epicenters outside domain that still contain enough presences (and absences)
if (neigh$epiType=='random') {

	points <- if (contrastSupplied) {
		SpatialPoints(coords=data.frame(rbind(pres[ , 1:2], contrast[ , 1:2])), proj4string=CRS(projection(domain)))
	} else {
		SpatialPoints(coords=data.frame(pres[ , 1:2]), proj4string=CRS(projection(domain)))
	}
	points <- spTransform(points, CRS(eaCRS))

	region <- extent(c(extent(points)@xmin - max(neigh$radius), extent(points)@xmax + max(neigh$radius), extent(points)@ymin - max(neigh$radius), extent(points)@ymax + max(neigh$radius)))
	
	region <- as(region, Class='SpatialPolygons')
	projection(region) <- eaCRS

}

##########
## MAIN ##
##########

### random placement of epicenters
##################################

if (neigh$epiType=='random') {

	# for each radius size
	for (radiusIndex in seq_along(neigh$radius)) {
	
		if (verbose) say('Placing neighborhoods of radius ', neigh$radius[radiusIndex], '.')
		
		# define region from which epicenter selection occurs (smaller for smaller radii)
		selRegion <- extent(c(extent(points)@xmin - neigh$radius[radiusIndex], extent(points)@xmax + neigh$radius[radiusIndex], extent(points)@ymin - neigh$radius[radiusIndex], extent(points)@ymax + neigh$radius[radiusIndex]))

		selRegion <- as(selRegion, Class='SpatialPolygons')
		projection(selRegion) <- eaCRS

		tries <- 0 # number of times trying to locate epicenters

		# progress bar
		if (verbose) pb <- txtProgressBar(min=0, max=neigh$numDesired[radiusIndex], width=30, style=3)
	
		# while still watning to try and number of neighborhoods is < number desired
		while (tries <= neigh$attempts & swarm$neigh$numNeigh[radiusIndex] < neigh$numDesired[radiusIndex]) {
		
			tries <- tries + 1
			
			# candidate epicenter site
			epiTry <- spsample(selRegion, n=neigh$numDesired[radiusIndex] - swarm$neigh$numNeigh[radiusIndex], type='random')
			epiTrySp <- sp::spTransform(epiTry, CRS(projection(domain)))
			
			# for each candidate neighborhood
			for (countNeigh in seq_along(epiTry)) {
			
				# see if contains enough presences
				distPres <- neigh$distFx(p1=epiTrySp[countNeigh], p2=presSp)
				presIn <- distPres <= neigh$radius[radiusIndex]
				numPresIn <- sum(presIn)
				weightPresIn <- sum(rawPresWeights[presIn])
				minPresMet <- numPresIn >= minPres
				minPresWeightMet <- weightPresIn >= minPresWeight
				
				# see if contains enough absences
				if (contrastSupplied) {
					distContrast <- neigh$distFx(p1=epiTrySp[countNeigh], p2=presSp)
					contrastIn <- distContrast <= neigh$radius[radiusIndex]
					numContrastIn <- sum(contrastIn)
					weightContrastIn <- sum(rawContrastWeights[contrastIn])
					minContrastMet <- numContrastIn >= minContrast
					minContrastWeightMet <- weightContrastIn >= minContrastWeight
				} else {
					numContrastIn <- NA
					minContrastMet <- TRUE
					minContrastWeightMet <- TRUE
				}
				
				# neighborhood has enough presences (and absences), remember it
				if (minPresMet & minContrastMet & minPresWeightMet & minContrastWeightMet) {

					setTxtProgressBar(pb, swarm$neigh$numNeigh[radiusIndex])

					# increment number of neighborhoods successfully selected
					swarm$neigh$numNeigh[radiusIndex] <- swarm$neigh$numNeigh[radiusIndex] + 1

					# calculate weight of each presence (and absence) and total neighborhood weight
					globalPresWeight <- neigh$rawNeighWeights[radiusIndex] * rawPresWeights[which(presIn)]
					if (neigh$weighting=='gauss') {
						globalPresWeight <- globalPresWeight * dnorm(distPres[presIn], mean=0, sd=neigh$sd[radiusIndex])
					}
					neighWeight <- sum(globalPresWeight)
					
					if (contrastSupplied) {
						globalContrastWeight <- neigh$rawNeighWeights[radiusIndex] * rawContrastWeights[which(contrastIn)]
						if (neigh$weighting=='gauss') globalContrastWeight <- globalContrastWeight * dnorm(distContrast[contrastIn], mean=0, sd=neigh$sd[radiusIndex])
						neighWeight <- neighWeight + sum(globalContrastWeight)
					} else {
						globalContrastWeight <- NA
					}

					# remember neighborhood stats
					swarm$neighStats <- rbind(
						swarm$neighStats,
						data.frame(
							type='delimited',
							weighting=neigh$weighting,
							radius=neigh$radius[radiusIndex],
							sd=if (exists('sd', where=neigh, inherits=FALSE)) { neigh$sd[radiusIndex] } else { NA },
							neighWeight=neighWeight,
							numPres=numPresIn,
							numContrast=numContrastIn,
							longitude=coordinates(epiTrySp[countNeigh])[1],
							latitude=coordinates(epiTrySp[countNeigh])[2]
						)
					)

					# remember index of presences (and absences) and weights in this neighborhood and select which are in/out of bag
					presIndex <- which(presIn)
					presIndexInBag <- sort(sample(presIndex, round(bag * length(presIndex))))
					presIndexOutOfBag <- presIndex[!(presIndex %in% presIndexInBag)]
					
					if (contrastSupplied) {
						contrastIndex <- which(contrastIn)
						contrastIndexInBag <- sort(sample(contrastIndex, round(bag * length(contrastIndex))))
						contrastIndexOutOfBag <- contrastIndex[!(contrastIndex %in% contrastIndexInBag)]
					} else {
						contrastIndex <- contrastIndexInBag <- contrastIndexOutOfBag <- NA
					}
					
					swarm$neighPresContrast[[length(swarm$neighPresContrast) + 1]] <- list(
						presIndex=presIndex,
						globalPresWeight=globalPresWeight,
						presIndexInBag=presIndexInBag,
						presIndexOutOfBag=presIndexOutOfBag,
						contrastIndex=contrastIndex,
						globalContrastWeight=globalContrastWeight,
						contrastIndexInBag=contrastIndexInBag,
						contrastIndexOutOfBag=contrastIndexOutOfBag
					)

					# remember coverage
					swarm$pres$coverage[presIndex] <- swarm$pres$coverage[presIndex] + 1
					if (contrastSupplied) swarm$contrast$coverage[contrastIndex] <- swarm$contrast$coverage[contrastIndex] + 1
					
				} # if sufficient presences (and abseces) in this neighborhood
				
			}  # next candidate neighborhood
			
		} # if still needing to place neighborhoods of this size

		if (verbose) say('Successfully placed ', swarm$neigh$numNeigh[radiusIndex], ' neighborhoods (', round(100 * swarm$neigh$numNeigh[radiusIndex] / neigh$numDesired[radiusIndex], 1), '% of desired).', pre=1)
		
	} # next radius
	
}

### omnibus models (cover entire area)
######################################

if (neigh$omnibus > 0) {

	# for each omnibus neighborhood
	for (countOmni in 1:neigh$omnibus) {

		# calculate weight of each presence (and absence) and total neighborhood weight
		neighWeight <- neigh$omnibusWeightScale * mean(swarm$neighStats$neighWeight[swarm$neighStats$type=='delimited'])
		globalPresWeight <- rawPresWeights
		
		globalContrastWeight <- if (contrastSupplied) {
			rawContrastWeights
		} else {
			NA
		}

		# remember neighborhood stats
		swarm$neighStats <- rbind(
			swarm$neighStats,
			data.frame(
				type='omnibus',
				weighting='user-supplied',
				radius=NA,
				sd=NA,
				neighWeight=neighWeight,
				numPres=nrow(pres),
				numContrast=if (contrastSupplied) { nrow(contrast) } else { NA },
				longitude=NA,
				latitude=NA
			)
		)

		# remember index of presences (and absences) in this neighborhood and in-bag/out-of-bag status (since this is an omnibus model all presences (and absences) are automatically in
		presIndex <- 1:nrow(pres)
		presIndexInBag <- sort(sample(presIndex, round(bag * length(presIndex))))
		presIndexOutOfBag <- presIndex[!(presIndex %in% presIndexInBag)]
		
		if (contrastSupplied) {
			contrastIndex <- 1:nrow(contrast)
			contrastIndexInBag <- sort(sample(contrastIndex, round(bag * length(contrastIndex))))
			contrastIndexOutOfBag <- contrastIndex[!(contrastIndex %in% contrastIndexInBag)]
		} else {
			contrastIndex <- contrastIndexInBag <- contrastIndexOutOfBag <- NA
		}
		
		swarm$neighPresContrast[[length(swarm$neighPresContrast) + 1]] <- list(
			presIndex=presIndex,
			globalPresWeight=globalPresWeight,
			presIndexInBag=presIndexInBag,
			presIndexOutOfBag=presIndexOutOfBag,
			contrastIndex=contrastIndex,
			globalContrastWeight=globalContrastWeight,
			contrastIndexInBag=contrastIndexInBag,
			contrastIndexOutOfBag=contrastIndexOutOfBag
		)

	} # next omnibus neighborhood

}

#####################
## POST-PROCESSING ##
#####################

if (any(swarm$pres$coverage == 0)) warning('Not all presences are covered by a delimited neighborhood!', immediate.=TRUE)
if (contrastSupplied) if (any(swarm$contrast$coverage == 0)) warning('Not all contrast sites are covered by a delimited neighborhood!', immediate.=TRUE)

return(swarm)
	
}

