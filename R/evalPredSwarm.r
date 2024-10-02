evalPredSwarm <- function(
	swarm,
	domain,
	modelDir,
	model=c('glm', 'gam', 'brt', 'randomForest'),
	metric=c('pred', 'auc', 'cbi', 'fpb', 'cor', 'msss_thold', 'msss_sens', 'msss_spec', 'mdss_thold', 'mdss_sens', 'mdss_spec', 'minTrainPred_thold', 'minTrainPred_sens', 'minTrainPred_spec', 'train10Perc_thold', 'train10Perc_sens', 'train10Perc_spec'),
	sampleSize=10000,
	incOmnibus=TRUE,
	verbose=TRUE,
	...
) {
# evalPredSwarm Produces evaluation raster stack for a swarm model.
#
# ARGUMENTS
# swarm			Object of class "swarm".  Produced by function makeSwarm().
# vars			Character list of variables to be used in model.
# modelDir		Character, directory in which to save models
# model			Character, model function to use for submodels.  Currently supports:
#					'brt' (gradient boosted machines/boosted regression trees)
#					'glm' (generalized linear model)
#					'gam' (generalized additive model)
# verbose		Logical, if TRUE then reports progress.
# ...			Arguments to pass to the model function including:
#					glm and gam: formula (e.g., "formula=as.formula(response ~ x + z + I(x^2)'"; note that the dependent variable should ALWAYS be named "response") and family (e.g., "family='binomial'")
#					gam: method, optimizer, scale, select, gamma, etc.
#					gbm.step: tree.complexity, learning.rate, bag.fraction, n.folds, prev.stratify, n.trees, step.size, max.trees, verbose, silent
# verbose		TRUE ==> Display detailed progress.  FALSE ==> Display minimal progress.
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
# SOURCE	source('C:/ecology/Drive/R/SDM/SSSDM/evalSwarm.r')
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

####################
## PRE-PROCESSING ##
####################

if (verbose) say('=== Evaluating Swarm Model ===')

# create output stack
template <- domain[[1]] * 0

evalStack <- stack(template)
if (length(metric) > 1) for (i in 2:length(metric)) evalStack <- stack(evalStack, template)
names(evalStack) <- metric

# add layers for weights and number of neighborhoods
evalStack <- stack(evalStack, template, template)
names(evalStack)[(nlayers(evalStack) - 1):nlayers(evalStack)] <- c('weight', 'numNeigh')

##########
## MAIN ##
##########

	### for each neighborhood ###
	#############################

	for (countNeigh in nrow(swarm$neighStats):1) {

		if (verbose) say('Evaluating/predicting swarm model #', countNeigh, ' of ', sum(swarm$neighStats$type=='delimited') + ifelse(incOmnibus, sum(swarm$neighStats$type=='omnibus'), 0), ' at ', date(), '.')

		# get model
		load(paste0(modelDir, '/', toupper(model), ' Submodel ', prefix(countNeigh, 4), '.Rdata'))

		# get test data
		presContrast <- submodel$allData$response
		trainTest <- submodel$allData$trainTest
		
		# predict
		preds <- predict(submodel$submodel, newdata=submodel$allData, type='response', n.trees=theModel$gbm.call$best.trees, na.rm=TRUE, ...)
		# preds <- predict(submodel$submodel, newdata=submodel$allData, type='response', n.trees=theModel$gbm.call$best.trees, na.rm=TRUE)
		
		predPresTrain <- preds[presContrast==1 & trainTest=='train']
		predContrastTrain <- preds[presContrast==0 & trainTest=='train']
		
		predPresTest <- preds[presContrast==1 & trainTest=='test']
		predContrastTest <- preds[presContrast==0 & trainTest=='test']
		
		# evaluations
		evalTrain <- evaluate(p=as.vector(predPresTrain), a=as.vector(predContrastTrain), tr=seq(0, 1, by=0.01))
		evalTest <- evaluate(p=as.vector(predPresTest), a=as.vector(predContrastTest), tr=seq(0, 1, by=0.01))
		
		### get (cropped) environmental stack and weights raster
		########################################################

		# delimited neighborhoods
		if (swarm$neighStats$type[countNeigh]=='delimited') {

			# get epicenter
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
		
		# omnibus models
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

		# remember neighborhood mask raster
		if (swarm$neighStats$type[countNeigh] == 'delimited') {
		
			targetIndex <- which(names(evalStack) %in% 'numNeigh')
			evalStack[[targetIndex]] <- mosaic(evalStack[[targetIndex]], domainWeightRast > 0, fun=sum)
			names(evalStack)[[targetIndex]] <- 'numNeigh'
		
		}

		# get random background sites and predictions for some metrics
		if (any(metric %in% c('auc', 'cbi', 'fpb', 'msss_thold', 'msss_spec', 'mdss_thold', 'mdss_spec', 'minTrainPred_thold', 'minTrainPred_spec', 'train10Perc_thold', 'train10Perc_spec'))) {
			
			randSites <- sampleRast(mask=domainNeigh, n=sampleSize, adjArea=TRUE, replace=TRUE, prob=FALSE)
			env <- as.data.frame(extract(domainNeigh, randSites))
			predBg <- predict(submodel$submodel, newdata=env, na.rm=TRUE, type='response',n.trees=submodel$submodel$gbm.call$best.trees)
			
			neighEval <- evaluate(p=as.vector(predPresTest), a=as.vector(predBg), tr=seq(0, 1, by=0.01))
			
		}
		
		### for each output metric
		##########################
		
		for (thisMetric in metric) {

			### prediction
			##############
			
			if (thisMetric == 'pred') {
				
				outRast <- if (model %in% c('glm', 'gam')) {
				
					predict(
						object=domainNeigh,
						model=submodel$submodel,
						na.rm=TRUE,
						type='response'
					)
					
				} else if (model=='brt') {
					
					predict(
						domainNeigh,
						submodel$submodel,
						n.trees=submodel$submodel$gbm.call$best.trees,
						na.rm=TRUE,
						type='response'
					)
					
				}
				
			}

			### other metrics
			#################
			
			if (thisMetric == 'weight') outRast <- domainTemplate + 1
			if (thisMetric == 'neighMask') outRast <- domainTemplate + 1
			
			if (thisMetric == 'auc') outRast <- domainTemplate + neighEval@auc
			if (thisMetric == 'cbi') outRast <- domainTemplate + contBoyce(predAtPres=predPresTest, predAtBg=predBg, numClasses=10)
			if (thisMetric == 'fpb') outRast <- domainTemplate + mean(calcFpb(predTestPres=predPresTest, predTestAbs=predBg, tr=seq(0, 1, by=0.01)), na.rm=TRUE)
			if (thisMetric == 'cor') outRast <- domainTemplate + neighEval@cor

			if (thisMetric == 'msss_thold') outRast <- domainTemplate + neighEval@t[which.max(neighEval@TPR + neighEval@TNR)]
			if (thisMetric == 'msss_sens') outRast <- domainTemplate + (sum(predPresTest >= neighEval@t[which.max(neighEval@TPR + neighEval@TNR)], na.rm=TRUE) / length(na.omit(predPresTest)))
			if (thisMetric == 'msss_spec') outRast <- domainTemplate + (sum(predBg < neighEval@t[which.max(neighEval@TPR + neighEval@TNR)], na.rm=TRUE) / length(na.omit(predBg)))
			
			if (thisMetric == 'mdss_thold') outRast <- domainTemplate + neighEval@t[which.min(abs(neighEval@TPR - neighEval@TNR))]
			if (thisMetric == 'mdss_sens') outRast <- domainTemplate + (sum(predPresTest >= neighEval@t[which.min(abs(neighEval@TPR - neighEval@TNR))], na.rm=TRUE) / length(na.omit(predPresTest)))
			if (thisMetric == 'mdss_spec') outRast <- domainTemplate + (sum(predBg < neighEval@t[which.min(abs(neighEval@TPR - neighEval@TNR))], na.rm=TRUE) / length(na.omit(predBg)))

			if (thisMetric == 'minTrainPred_thold') outRast <- domainTemplate + min(predPresTrain, na.rm=TRUE)
			if (thisMetric == 'minTrainPred_sens') outRast <- domainTemplate + sum(predPresTest >= min(predPresTrain, na.rm=TRUE), na.rm=TRUE) / length(na.omit(predPresTest))
			if (thisMetric == 'minTrainPred_spec') outRast <- domainTemplate + sum(predBg < min(predPresTrain, na.rm=TRUE), na.rm=TRUE) / length(na.omit(predBg))
			
			if (thisMetric == 'train10Perc_thold') outRast <- domainTemplate + quantile(predPresTrain, 0.1, na.rm=TRUE)
			if (thisMetric == 'train10Perc_sens') outRast <- domainTemplate + sum(predPresTest >= quantile(predPresTrain, 0.1, na.rm=TRUE), na.rm=TRUE) / length(na.omit(predPresTest))
			if (thisMetric == 'train10Perc_spec') outRast <- domainTemplate + sum(predBg < quantile(predPresTrain, 0.1, na.rm=TRUE), na.rm=TRUE) / length(na.omit(predBg))
			
			### weigh
			#########
			
			outRast <- outRast * domainWeightRast
			
			### insert results into stack
			#############################
			
			targetIndex <- which(names(evalStack) %in% thisMetric)
			evalStack[[targetIndex]] <- mosaic(evalStack[[targetIndex]], outRast, fun=sum)
			names(evalStack)[[targetIndex]] <- thisMetric
			
			rm(outRast)
			
		} # next evaluation metric
		
		rm(domainTemplate, domainWeightRast, domainNeigh)
				
	} # next neighborhood

#####################
## POST-PROCESSING ##
#####################

	### generate neighborhood mask
	numNeighIndex <- which(names(evalStack) %in% 'numNeigh')
	numNeigh <- evalStack[[numNeighIndex]]
	neighMask <- calc(numNeigh, fun=function(x, ...) ifelse(x > .Machine$double.eps, 1, NA), na.rm=TRUE)
	evalStack <- stack(evalStack, neighMask)
	names(evalStack)[nlayers(evalStack)] <- 'neighMask'

	### calculate weighted average for each evaluation metric/prediction
	toWeighIndex <- {1:length(metric)}[!(metric %in% c('weight', 'neighMask', 'numNeigh'))]
	weightRast <- evalStack[[which(names(evalStack) %in% 'weight')]]
	for (i in toWeighIndex) evalStack[[i]] <- evalStack[[i]] / weightRast
	names(evalStack)[1:length(metric)] <- metric

	if (verbose) say('DONE!', pre=1)

	return(evalStack)

}

