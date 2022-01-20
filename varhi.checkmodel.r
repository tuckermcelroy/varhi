varhi.checkmodel <- function(data,phi.arry,var.pred,N.final,shapiro=TRUE,subbarao=TRUE,H=0)
{

	########################################################################
	#	varhi.checkmodel by Tucker McElroy, October 17, 2017
	#
	#    Purpose: given a VAR model with final core and auxiliary variables,
	#	fit and examine residuals for goodness of fit.
 	#    Inputs: 
	#	data: T x N  matrix of N  time series of length T
	#	phi.arry: array of N x N x p.mdl of the fitted VAR coefficients
	#	var.pred: N x N covariance matrix for innovations
	#	N.final: indices of final variables
	#	shapiro: flag for doing the shapiro-wilk normality test
	#	subbarao: flag for doing the subbarao stationarity test
	#	H: number of lags for portmanteau, set to zero if not desired
	#    Outputs: print to screen up to three diagnostics
	#	returns array of N x N x p.mdl of standard errors for fitted VAR
	#    Requires: package "fractal"
	#
	########################################################################		

data <- data[,N.final] 
T <- dim(data)[1]
N <- dim(data)[2]
mu <- colMeans(data)
data <- t(t(data) - mu)

p.mdl <- dim(phi.arry)[[3]]
phi.mat <- matrix(phi.arry,nrow=N,ncol=N*p.mdl)
data.mat <- t(data)[,(p.mdl+1):T]
for(i in 1:p.mdl) { data.mat <- rbind(data.mat,t(data)[,(p.mdl+1-i):(T-i)]) }
resids <- cbind(diag(N),-phi.mat) %*% data.mat
resids <- t(resids)
gamma.acf <- acf(resids,na.action=na.pass,lag.max=H,plot=FALSE)

####################
### Model checking
 
# normality
if(shapiro) {
for(i in 1:N) { print(shapiro.test(resids[,i])) } }
 
# stationarity
if(subbarao) {
for(i in 1:N) { print(summary(stationarity(resids[,i]))$test) } }

# Portmanteau
if(H > 0) {
all.port <- 0
ports <- rep(0,N)
for(i in 1:H) 
{
	all.port <- all.port + sum(diag(t(gamma.acf$acf[i+1,,]) %*% gamma.acf$acf[i+1,,] ))
	ports <- ports + diag(gamma.acf$acf[i+1,,])^2
}
all.port <- T*all.port
ports <- T*ports
print(1-pchisq(all.port,df= N^2*(H-p.mdl)))
print(1-pchisq(ports,df= H-p.mdl))
}

#####################
## standard errors

mdl.acf <- VARMAauto(NULL,phi.arry,var.pred,p.mdl)

## get Toeplitz covariance matrices
if(p.mdl==1) 
{
	mdl.toep <- mdl.acf[,,1]
}
if(p.mdl==2) 
{
	mdl.toep <- rbind(cbind(mdl.acf[,,1],mdl.acf[,,2]),
		cbind(t(mdl.acf[,,2]),mdl.acf[,,1]))
}
if(p.mdl==3) 
{
	mdl.toep <- rbind(cbind(mdl.acf[,,1],mdl.acf[,,2],mdl.acf[,,3]),
			cbind(t(mdl.acf[,,2]),mdl.acf[,,1],mdl.acf[,,2]),
			cbind(t(mdl.acf[,,3]),t(mdl.acf[,,2]),mdl.acf[,,1]))
}
if(p.mdl==4) 
{
	mdl.toep <- rbind(cbind(mdl.acf[,,1],mdl.acf[,,2],mdl.acf[,,3],mdl.acf[,,4]),
			cbind(t(mdl.acf[,,2]),mdl.acf[,,1],mdl.acf[,,2],mdl.acf[,,3]),
			cbind(t(mdl.acf[,,3]),t(mdl.acf[,,2]),mdl.acf[,,1],mdl.acf[,,2]),
			cbind(t(mdl.acf[,,4]),t(mdl.acf[,,3]),t(mdl.acf[,,2]),mdl.acf[,,1])) 
}

asymp.var <- solve(mdl.toep) %x% var.pred
asymp.se <- sqrt(diag(asymp.var))/sqrt(T)
asymp.se <- array(matrix(asymp.se,nrow=N),c(N,N,p.mdl))

return(list(asymp.se,resids))

}