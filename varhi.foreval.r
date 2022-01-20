varhi.foreval <- function(data,phi.arry,N.final)
{

	########################################################################
	#	varhi.foreval by Tucker McElroy, August 15, 2017
	#
	#    Purpose: given a VAR model with final core and auxiliary variables,
	#	examine forecast performance
 	#    Inputs: 
	#	data: T x N  matrix of N time series of length T
	#	phi.arry: array of N x N x p.mdl of the fitted VAR coefficients
	#	N.final: indices of final variables
	#    Outputs: data.fore is T x N.final matrix of 1-step ahead forecasts,
	#	for each of the N.final series, based on the best VAR model
	#
	########################################################################		

T.total <- dim(data)[1]
data.final <- data[,N.final]
mu <- colMeans(data.final)
data.final <- t(t(data.final) - mu)
p.mdl <- dim(phi.arry)[3]

if(p.mdl==1) { 
	phi.hat <- phi.arry[,,1]
}
if(p.mdl==2) {
	phi1.hat <- phi.arry[,,1]
	phi2.hat <- phi.arry[,,2]
}
if(p.mdl==3) {
	phi1.hat <- phi.arry[,,1]
	phi2.hat <- phi.arry[,,2]
	phi3.hat <- phi.arry[,,3]
}
if(p.mdl==4) {
	phi1.hat <- phi.arry[,,1]
	phi2.hat <- phi.arry[,,2]
	phi3.hat <- phi.arry[,,3]
	phi4.hat <- phi.arry[,,4]
}

start.eval <- p.mdl
fore.lead <- T.total-start.eval
if(fore.lead > 0) {
forecasts <- NULL
for(h in 1:fore.lead)
{
	if(p.mdl==1) { forecast <- mu + phi.hat %*% (data[(start.eval+h-1),N.final] - mu) }
	if(p.mdl==2) { forecast <- mu + phi1.hat %*% (data[(start.eval+h-1),N.final] - mu) +
				phi2.hat %*% (data[(start.eval+h-2),N.final] - mu) }
	if(p.mdl==3) { forecast <- mu + phi1.hat %*% (data[(start.eval+h-1),N.final] - mu) +
				phi2.hat %*% (data[(start.eval+h-2),N.final] - mu) +
				phi3.hat %*% (data[(start.eval+h-3),N.final] - mu)}
	if(p.mdl==4) { forecast <- mu + phi1.hat %*% (data[(start.eval+h-1),N.final] - mu) +
				phi2.hat %*% (data[(start.eval+h-2),N.final] - mu) +
				phi3.hat %*% (data[(start.eval+h-3),N.final] - mu) +
				phi4.hat %*% (data[(start.eval+h-4),N.final] - mu) }
	forecasts <- cbind(forecasts,forecast)
} 
data.fore <- rbind(data[1:start.eval,N.final],t(forecasts)) 
}	# end if

return(data.fore)
}

