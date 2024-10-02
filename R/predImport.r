predImport <- function(
	swarm,
	domain,
	modelDir,
	model=c('glm', 'gam', 'brt', 'randomForest'),
	metric=c('importRPresContrast', 'importRBg', 'importRBgStratScaled', 'importCbi', 'importCbiScaled', 'importAuc', 'importCor'),
	incOmnibus=TRUE,
	perms=100,
	sampleSize=10000,
	verbose=TRUE,
	...
) {
# predImport Evaluates importance of predictors in a swarm model.
#
# ARGUMENTS
# swarm			Object of class "swarm".  Produced by function makeSwarm().
# vars			Character list of variables to be used in model.
# modelDir		Character, directory in which to save models
# model			Character, model function to use for submodels.  Currently supports:
#					'brt' (gradient boosted machines/boosted regression trees)
#					'glm' (generalized linear model)
#					'gam' (generalized additive model)
# metric		Character list, indcates what type(s) of importance metric to be calculated. The output will have one raster per metric per variable (e.g., if two metrics are used then the output will have 2 * <number of variables> rasters).  The "correlation" metric used below is the correlation between observed predictions at test sites with predictions at the same sites with each variable permuted in turn.  For all metrics lower values connote greater importance. Also, for each metric values are calcualted perm times per neighborhood for each variable then averaged.
#					importRPresContrast		Correlation metric with test sites being test presences and test contrast sites.
#					importRBg	Correlation metric with test sites being randomly located background sites (equal in number to to sampleSize).
#					importRBgStrat Correlation metric with test sites being background sites drawn in a stratified manner from the prediction raster. The function attempts to calculate 10 bins of predicted values (from 0 to 0.1, 0.1 to 0.2, 0.2 to 0.3, ..., 0.9 to 1). It then samples from each bin round(sampleSize  / 10) sites at which predictions are made (with environmental values as observed or permuted).
#					importRBgStratScaled Correlation metric with test sites being background sites drawn in a stratified manner from the prediction raster *after the prediction raster has been rescaled to [0, 1]. The function calculates 10 bins of rescaled predicted values (from 0 to 0.1, 0.1 to 0.2, 0.2 to 0.3, ..., 0.9 to 1). It then samples from each bin round(sampleSize  / 10) sites at which predictions are made (with environmenta values as observed or permuted). In contrast to importRBgStrat, this scheme is more senistive to effects of predictors on gradients in suitability, even if the gradients are shallow.
#					importCbi	Continuous Boyce Index calculated using test presences and sampleSize randomly located background sites. 10 bins are used spanning from the minimum to maxim predicted value in the neighborhood.  In contrast to importCbiScaled this metric emphasizes importance of predictors across gradients of suitability, even if they are shallow.
#					importCbiScaled Continous Boyce Index calculated using sampleSize randomly located background sites.  10 bins are used spanning the range [0, 1]. In contrast to importCbi, this metric reflects absolute importance.
#					importAuc	AUC calculated using test presences and sampleSize randomly located background sites.??????????
#					importCor	COR calculated using test presences and sampleSize randomly located background sites.??????????
# incOmnibus	Logical, if TRUE then omnibus models are included in importance calculations.  If FALSE then only delimited neighborhood models are used.
# perms			Positive integer, number of times to permute variable of interest then calculate response metric. For each neighborhood and metric the average is taken over these values.
# sampleSize	Positive integer, number of randomly located background sites to use for calculation of metrics.
# verbose		Logical, if TRUE then reports progress.
# ...			Additional arguments (not implemented).
#
# VALUES
# 
# 
# REQUIRED DEPENDANCIES
#
# OPTIONAL DEPENDANCIES
# dismo, mgcv
#
# BAUHAUS
# 
#
# EXAMPLE
# FUNCTION()
#
# SOURCE	source('C:/ecology/Drive/R/SDM/SSSDM/predImport.r')
#
# TESTING
#
#
# LICENSE
# This document is copyright ©2016 by Adam B. Smith.  This document is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3, or (at your option) any later version.  This document is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. Copies of the GNU General Public License versions are available at http://www.R-project.org/Licenses/.
#
# AUTHOR	Adam B. Smith | Missouri Botanical Garden, St. Louis, Missouri | adamDOTsmithATmobotDOTorg
# DATE		2016-12
# REVISIONS 

############################
## FUNCTIONS AND PACKAGES ##
############################

model <- match.arg(model)

if (model=='gam') library(mgcv)
if (model=='brt') library(dismo)

if ('importCbi' %in% metric) source('C:/ecology/Drive/R/SDM/SSSDM/contBoyce.r')

####################
## PRE-PROCESSING ##
####################

	if (verbose) say('=== Calculating Predictor Importance of Swarm Model ===')

	### create output stack
	#######################
	
	template <- domain[[1]] * 0

	evalStack <- stack(template)
	for (i in 2:nlayers(domain)) evalStack <- stack(evalStack, template)
	if (length(metric) > 1) for (i in 1:((length(metric) - 1) * nlayers(domain))) evalStack <- stack(evalStack, template)

	names(evalStack) <- paste0(rep(metric, each=nlayers(domain)), '_', names(domain))

	evalStack <- stack(evalStack, template)
	names(evalStack)[nlayers(evalStack)] <- 'weight'

##########
## MAIN ##
##########

	#############################
	### for each neighborhood ###
	#############################

	if (verbose) pb <- txtProgressBar(min=0, max=nrow(swarm$neighStats), width=30, style=3)
		
	for (countNeigh in 1:nrow(swarm$neighStats)) {

		# if (verbose)  say('Evaluating predictor importance of swarm model #', countNeigh, ' of ', sum(swarm$neighStats$type=='delimited') + ifelse(incOmnibus, sum(swarm$neighStats$type=='omnibus'), 0), ' at ', date(), '.', post=0)
		# get model
		load(paste0(modelDir, '/', toupper(model), ' Submodel ', prefix(countNeigh, 4), '.Rdata'))

		# get environmental data frames for test presences/contrasts
		envTestPres <- submodel$allData[submodel$allData$response == 1 & submodel$allData$trainTest == 'train', ]
		envTestContrast <- submodel$allData[submodel$allData$response == 0 & submodel$allData$trainTest == 'train', ]
		
		# predict... add small random amount to avoid errors from sd(x) = 0
		predTestPres <- predict(submodel$submodel, newdata=envTestPres, type='response', n.trees=theModel$gbm.call$best.trees, na.rm=TRUE) + rnorm(nrow(envTestPres), eps, eps)
		predTestContrast <- predict(submodel$submodel, newdata=envTestContrast, type='response', n.trees=theModel$gbm.call$best.trees, na.rm=TRUE) + rnorm(nrow(envTestContrast), eps, eps)
		
		# evaluations
		evalTest <- evaluate(p=as.vector(predTestPres), a=as.vector(predTestContrast), tr=seq(0, 1, by=0.01))
		
		### get (cropped) environmental stack and weights raster
		########################################################
		
		# delimited neighborhoods
		if (swarm$neighStats$type[countNeigh]=='delimited') {

			# epicenter
			epi <- swarm$neighStats[countNeigh, c('longitude', 'latitude')]
		
			# crop/mask domain to neighborhood
			domainNeigh <- cropToNeighborhood(x=domain, epi=epi, radius=swarm$neighStats$radius[countNeigh], CRS=swarm$CRS, eaCRS=swarm$eaCRS, quadsegs=20)
			
			# template to which to add/subtract performance values, weights, etc.
			domainTemplate <- domainNeigh[[1]] * 0
			names(domainTemplate) <- 'domainTemplate'

			# if weighting by distance to neighborhood epicenter
			if (swarm$neighStats$weighting[countNeigh]=='gauss') gaussWeight <- weighByGauss(x=domainTemplate, epi=epi, sd=swarm$neighStats$sd[countNeigh], distFx=swarm$neigh$distFx)

			# weights raster
			domainWeightRast <- domainTemplate + swarm$neighStats$neighWeight[countNeigh]
			if (swarm$neighStats$weighting[countNeigh]=='gauss') domainWeightRast <- domainWeightRast * gaussWeight

		## omnibus models
		} else if (swarm$neighStats$type[countNeigh]=='omnibus') {
		
			domainNeigh <- domain

			# template to which to add/subtract performance values, weights, etc.
			domainTemplate <- domainNeigh[[1]] * 0
			names(domainTemplate) <- 'domainTemplate'

			domainWeightRast <- domainTemplate + swarm$neighStats$neighWeight[countNeigh]
			
		}

		# remember weights raster
		targetIndex <- which(names(evalStack) %in% 'weight')
		evalStack[[targetIndex]] <- mosaic(evalStack[[targetIndex]], domainWeightRast, fun=sum)
		names(evalStack)[[targetIndex]] <- 'weight'

		### prediction
		##############

		if (any(metric %in% c('importRBgStrat', 'importRBgStratScaled'))) {
			
			predRast <- predict(domainNeigh, submodel$submodel, n.trees=submodel$submodel$gbm.call$best.trees, na.rm=TRUE, type='response')
			
		}

		### for each output metric
		##########################
		
		for (thisMetric in metric) {

			# say(thisMetric, post=ifelse(thisMetric == metric[length(metric)], 1, 0))
			
			### generate data frame of predictors to be permuted... source of data depends on kind of test metric
			if (thisMetric %in% c('importRBg', 'importCbi', 'importAuc', 'importCor')) {
				
				# just randomly located sites
				testSites <- sampleRast(domainNeigh, n=sampleSize, adjArea=TRUE, replace=TRUE, prob=FALSE)
				
			} else if (thisMetric == 'importRBgStrat') {
			
				# NOT controlling for distribution of Pr(occ)
				predRastInteger <- cut(predRast, breaks=c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1))
				testSites <- sampleRastStrat(predRastInteger, nEach=round(sampleSize / 10), adjArea=TRUE, replace=TRUE)
				
			} else if (thisMetric == 'importRBgStratScaled') {

				# stratified random sample of sites after stretching raster
			
				# controlling for distribution of predicetd Pr(occ)
				predRastInteger <- stretch(predRast, 0, 1)
				predRastInteger <- cut(predRast, breaks=c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1))
				
				testSites <- sampleRastStrat(predRastInteger, nEach=round(sampleSize / 10), adjArea=TRUE, replace=TRUE)

			} else if (thisMetric == 'importRPresContrast') {
			
				# test preseneces and test contrast sites
				testSites <- rbind(envTestPres[ , 1:2], envTestContrast[ , 1:2])
				
			}
				
			# environmental data frame
			env <- data.frame(raster::extract(domainNeigh, testSites))
			env <- env[ , sort(names(env))]
			
			siteWeights <- if (swarm$neighStats$weighting[countNeigh]=='gauss') {
				c(raster::extract(gaussWeight, testSites))
			} else {
				rep(1, nrow(testSites))
			}
		
			### for each variable, calculate permutation importance
			#######################################################
			
			for (thisVar in names(domain)) {
				
				# assign data frame to contain permuted environmental values
				envPerm <- env

				# stores values of test statistic for each permutation
				out <- rep(NA, perms)
					
				# correlation permutation importance
				if (thisMetric %in% c('importRPresContrast', 'importRBg', 'importRBgStrat', 'importRBgStratScaled')) {
				
					predObs <- predict(submodel$submodel, newdata=env, type='response', n.trees=theModel$gbm.call$best.trees, na.rm=TRUE)
					
					# for each permutation
					for (i in 1:perms) {
					
						# permuted version of environmental data frame
						envPerm[ , thisVar] <- env[sample(1:nrow(env), nrow(env)), thisVar]
						predPerm <- predict(submodel$submodel, newdata=envPerm, type='response', n.trees=theModel$gbm.call$best.trees, na.rm=TRUE)
						out[i] <- cor(logitAdj(siteWeights * predObs, eps), logitAdj(siteWeights * predPerm, eps), use='na.or.complete')
					
					}
					
				} else if (thisMetric %in% c('importCbi', 'importAuc', 'importCor')) {
				
					# combine environmental data for presences/background
					thisPerm <- rbind(envTestPres[ , sort(names(domain))], envPerm)
					
					# for each permutation
					for (i in 1:perms) {
					
						# permuted version of environmental data frame
						thisPerm[ , thisVar] <- thisPerm[sample(1:nrow(thisPerm), nrow(thisPerm)), thisVar]
						
						# predict
						predPerm <- predict(submodel$submodel, newdata=thisPerm, type='response', n.trees=theModel$gbm.call$best.trees, na.rm=TRUE)
						
						predAtPres <- predPerm[1:nrow(envTestPres)]
						predAtBg <- predPerm[(nrow(envTestPres) + 1):length(predPerm)]
						
						# evaluation metric
						out[i] <- if (thisMetric == 'importCbi') {
							contBoyce(predAtPres=predAtPres, predAtBg=predAtBg, numClasses=10, presWeight=siteWeights[seq_along(predAtPres)], bgWeight=siteWeights[-seq_along(predAtPres)])
						} else if (thisMetric == 'importCbiScaled') {
							contBoyce(predAtPres=predAtPres, predAtBg=predAtBg, numClasses=10, presWeight=siteWeights[seq_along(predAtPres)], bgWeight=siteWeights[-seq_along(predAtPres)], autoWindow=FALSE)
						} else if (thisMetric == 'importAuc') {
							evaluate(p=as.vector(predAtPres), a=as.vector(predAtBg), tr=seq(0, 1, by=0.01))@auc
							aucWeighted(p=predAtPres, a=predAtBg, pw=siteWeights[seq_along(predAtPres)], aw=siteWeights[-seq_along(predAtPres)], na.rm=TRUE)
						} else if (thisMetric == 'importCor') {
							cor(siteWeights * c(predAtPres, predAtBg), c(rep(1, length(predAtPres)), rep(0, length(predAtBg))), use='na.or.complete')
							# evaluate(p=as.vector(predAtPres), a=as.vector(predAtBg), tr=seq(0, 1, by=0.01))@cor
						}
					
					}
				
				}
				
				### generate test metric raster
				out <- mean(out, na.rm=TRUE)
				if (is.na(out)) out <- 0
				
				outRast <- domainTemplate + out
				outRast <- outRast * domainWeightRast
			
				### insert results into stack
				#############################
				
				targetIndex <- which(names(evalStack) %in% paste0(thisMetric, '_', thisVar))
				evalStack[[targetIndex]] <- mosaic(evalStack[[targetIndex]], outRast, fun=sum)
				names(evalStack)[[targetIndex]] <- paste0(thisMetric, '_', thisVar)
				
				rm(outRast)
			
			} # next variable
			
		} # next metric
		
		rm(domainTemplate, domainWeightRast, domainNeigh)
		
		setTxtProgressBar(pb, countNeigh - 1)
		
	} # next neighborhood

#####################
## POST-PROCESSING ##
#####################

### number of neighborhoods raster and convert neighborhood mask to 1/NA
########################################################################

# calculate weighted average for each evaluation metric/prediction
toWeighIndex <- {1:(length(metric) * nlayers(domain))}[!(metric %in% c('weight'))]
weightRast <- evalStack[[which(names(evalStack) %in% 'weight')]]

for (i in toWeighIndex) {
	thisName <- names(evalStack[[i]])
	evalStack[[i]] <- evalStack[[i]] / weightRast
	names(evalStack)[i] <- thisName
}

if (verbose) {
	setTxtProgressBar(pb, countNeigh)
	say('')
}

evalStack

}

