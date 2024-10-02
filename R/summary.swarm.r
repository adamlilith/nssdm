summary.swarm <- function(swarm) {

	len <- 50 # lengh of a line before stats reported

	say('swarm object:', pre=1)
	say(suffix('names in domain', len=len), paste(swarm$domainNames, collapse=' '))
	say(suffix('CRS', len=len), swarm$CRS)
	say(suffix('equal-area CRS (eaCRS)', len=len), swarm$eaCRS)
	say(suffix('number of presences', len=len), nrow(swarm$pres))
	say(suffix('number of contrast sites', len=len), ifelse(is.null(swarm$contrast), NA, nrow(swarm$contrast)))
	say(suffix('minimum presences / neighborhood (user-defined)', len=len), swarm$minPres)
	say(suffix('minimum contrasts / neighborhood (user-defined)', len=len), swarm$minContrast)
	say(suffix('median presences / neighborhood', len=len), round(median(swarm$neighStats$numPres), 2))
	say(suffix('median contrasts / neighborhood', len=len), round(median(swarm$neighStats$numContrast), 2))
	say(suffix('presences outside delimited neighborhoods', len=len), sum(swarm$pres$coverage == 0))
	say(suffix('min/median/max presence coverage', len=len), paste(round(min(swarm$pres$coverage), 3), round(median(swarm$pres$coverage), 3), round(max(swarm$pres$coverage), 3), sep=' / '))
	if (!is.null(swarm$contrast)) {
		say(suffix('contrasts outside delimited neighborhoods', len=len), sum(swarm$contrast$coverage == 0))
		say(suffix('min/median/max contrast coverage', len=len), paste(round(min(swarm$contrast$coverage), 3), round(median(swarm$contrast$coverage), 3), round(max(swarm$contrast$coverage), 3), sep=' / '))
	}
	say(suffix('minimum total presence weight / neighborhood', len=len), swarm$minPresWeight)
	say(suffix('minimum total contrast weight / neighborhood', len=len), swarm$minContrastWeight)
	say(suffix('first few raw presence weights', len=len), paste(round(head(swarm$rawPresWeights), 2), collapse=' '))
	say(suffix('first few raw contrast weights', len=len), paste(round(head(swarm$rawContrastWeights), 2), collapse=' '))
	say(suffix('bag', len=len), swarm$bag)
	say(suffix('neighborhood radius', len=len), paste(swarm$neigh$radius, collapse=' '))
	if (swarm$neigh$weighting=='gauss') say(suffix('neighborhood sd', len=len), paste(swarm$neigh$sd, collapse=' '))
	say(suffix('raw neighborhood weighting', len=len), swarm$neigh$weighting)
	say(suffix('total neighborhoods', len=len), sum(swarm$neigh$numNeigh))
	say(suffix('number of actual neighborhoods', len=len), paste(swarm$neigh$numNeigh, collapse=' '))
	say(suffix('number of desired neighborhoods', len=len), paste(swarm$neigh$numDesired, collapse=' '))
	say(suffix('number of actual / desired neighborhoods', len=len), paste(round(swarm$neigh$numNeigh / swarm$neigh$numDesired, 2), collapse=' '))
	say(suffix('attempts', len=len), swarm$neigh$attempts)
	say(suffix('number of omnibus neighborhoods', len=len), sum(swarm$neighStats$type=='omnibus'))
	say(suffix('omnibus weighting scale', len=len), swarm$neigh$omnibusWeightScale)
	say('')
	
}
