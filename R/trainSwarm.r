trainSwarm <- function(
	swarm,
	vars,
	modelDir,
	model=c('glm', 'gam', 'brt', 'randomForest'),
	modelCon=TRUE,
	modelSel=TRUE,
	verbose=TRUE,
	...
) UseMethod('trainSwarm')


trainSwarm.default <- function(
	swarm,
	vars,
	modelDir,
	model=c('glm', 'gam', 'brt', 'randomForest'),
	modelCon=TRUE,
	modelSel=TRUE,
	verbose=TRUE,
	...
) {
# trainSwarm Trains a series of spatially distributed species distribution models.  Each model is trained using presence and contrast sites from spatial subsets of the entire study region.  These models can then be combined to produce non-stationary predictions across a geographic domain.
#
# ARGUMENTS
# swarm			Object of class "swarm".  Produced by function makeSwarm().
# vars			Character list of variables to be used in model.
# modelDir		Character, directory in which to save models
# model			Character, model function to use for submodels.  Currently supports:
#					'brt' (gradient boosted machines/boosted regression trees)
#					'glm' (generalized linear model)
#					'gam' (generalized additive model)
# modelSel		Logical. For some model algorithms model construction can be performed:
#					glm, gam (see ...)
# modelSel		Logical. For some model algorithms model selection can be performed:
#					glm, gam, maxent: AICc-based model selection (see ...)
# verbose		Logical, if TRUE then reports progress.
# ...			Arguments to pass to the model function including:
#					glm: see ?glm or for model construction/selection ?trainGlm
#					gam: see ?gam or for model construction/selection ?trainGam
#					gbm.step: see ?gbm.step
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
# SOURCE	source('C:/ecology/Drive/R/SDM/SSSDM/trainSwarm.r')
#
# TESTING
#
#
# LICENSE
# This document is copyright ©2016 by Adam B. Smith.  This document is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3, or (at your option) any later version.  This document is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. Copies of the GNU General Public License versions are available at http://www.R-project.org/Licenses/.
#
# AUTHOR	Adam B. Smith | Missouri Botanical Garden, St. Louis, Missouri | adamDOTsmithATmobotDOTorg
# DATE		2015-03-24
# REVISIONS 

############################
## FUNCTIONS AND PACKAGES ##
############################

model <- match.arg(model)

if (model=='gam') library(mgcv)
if (model=='brt') library(dismo)
if (modelSel) library(MuMIn)

####################
## PRE-PROCESSING ##
####################

if (verbose) say('=== Training Swarm Model ===')
dirCreate(modelDir)

ellipses <- list(...)

##########
## MAIN ##
##########

##########################
### for each epicenter ###
##########################
for (countNeigh in 1:nrow(swarm$neighStats)) {

	if (verbose) {
		# pb <- txtProgressBar(min=0, max=nrow(swarm$neighStats), width=30, style=3)
		say('Training model for neighborhood ', countNeigh, ' of ', nrow(swarm$neighStats), ' starting ', date(), '.')
	}

	# extract variables from swarm for ease of use
	presIndex <- swarm$neighPresContrast[[countNeigh]]$presIndex
	contrastIndex <- swarm$neighPresContrast[[countNeigh]]$contrastIndex

	presIndexInBag <- swarm$neighPresContrast[[countNeigh]]$presIndexInBag
	contrastIndexInBag <- swarm$neighPresContrast[[countNeigh]]$contrastIndexInBag

	presIndexOutOfBag <- swarm$neighPresContrast[[countNeigh]]$presIndexOutOfBag
	contrastIndexOutOfBag <- swarm$neighPresContrast[[countNeigh]]$contrastIndexOutOfBag

	pres <- swarm$pres[presIndex, ]
	contrast <- swarm$contrast[contrastIndex, ]

	pres$globalWeight <- swarm$neighPresContrast[[countNeigh]]$globalPresWeight
	contrast$globalWeight <- swarm$neighPresContrast[[countNeigh]]$globalContrastWeight

	pres$response <- 1
	contrast$response <- 0

	pres$trainTest <- ifelse(presIndex %in% presIndexInBag, 'train', 'test')
	contrast$trainTest <- ifelse(contrastIndex %in% contrastIndexInBag, 'train', 'test')

	allData <- rbind(pres, contrast)

	# (adding dummy variable for BRTs)
	if (model=='brt' & length(vars) == 1) {
		allData$DUMMY <- 1
		vars <- c(vars, 'DUMMY')
	}
	
	### manage weights
	##################

	### want total presence weight = total contrast weight and that no single weight is > 1
	
	modelPresWeight <- pres$globalWeight
	modelContrastWeight <- contrast$globalWeight
	
	sumPresWeight <- sum(modelPresWeight)
	sumContrastWeight <- sum(modelContrastWeight)

	# ensure equal weighting between presences and contrast sites
	if (sumPresWeight > sumContrastWeight) {
		modelContrastWeight <- modelContrastWeight * (sumPresWeight / sumContrastWeight)
	} else if (sumPresWeight < sumContrastWeight) {
		modelPresWeight <- modelPresWeight * (sumContrastWeight / sumPresWeight)
	}
	
	# ensure no single weight is > 1
	maxWeight <- max(modelPresWeight, modelContrastWeight)
	
	if (maxWeight > 1) {
		modelPresWeight <- modelPresWeight / maxWeight
		modelContrastWeight <- modelContrastWeight / maxWeight
	}

	allData$modelWeights <- c(modelPresWeight, modelContrastWeight)
	
	trainData <- allData[allData$trainTest == 'train', ]
	
	siteWeights <<- trainData$modelWeights / max(trainData$modelWeights) # declaring to global to avoid error in dredge() for glm, gam
	
	### train model
	###############
	
	# remember metadata
	submodel <- list()
	submodel$neighStats <- swarm$neighStats[countNeigh, ]
	submodel$vars <- vars
	submodel$allData <- allData
	
	# train model
	if (model=='brt') {
	
		submodel$submodel <- trainBrt(
			data=trainData,
			family='bernoulli',
			silent=TRUE,
			verbose=TRUE,
			w=siteWeights,
			...
		)
				
	} else if (model == 'gam') {
	
		submodel$submodel <- trainGam(
			resp='response',
			preds=vars,
			data=trainData,
			family='binomial',
			modelCon=modelCon,
			modelSel=modelSel,
			w=siteWeights,
			verbose=FALSE,
			...
		)
	
	} else if (model == 'glm') {

		submodel$submodel <- trainGlm(
			resp='response',
			preds=vars,
			data=trainData,
			modelCon=modelCon,
			modelSel=modelSel,
			family='binomial',
			w=siteWeights,
			verbose=FALSE,
			...
		)
	}
	
	### model selection
	###################

	# save
	save(submodel, file=paste0(modelDir, '/', toupper(model), ' Submodel ', prefix(countNeigh, 4), '.Rdata'), compress=TRUE)
	
	Sys.sleep(0.5)
	rm(submodel); gc()
	
	# setTxtProgressBar(pb, countNeigh)
	
} # next neighborhood

#####################
## POST-PROCESSING ##
#####################

if (verbose) say('DONE!')

}
