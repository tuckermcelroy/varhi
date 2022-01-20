varhi.fitmodel <- function(data,p.mdl,N.final,alpha,thresh)
{

	########################################################################
	#	varhi.fitmodel by Tucker McElroy, October 17, 2017
	#
	#    Purpose: given a VAR model with final core and auxiliary variables,
	#	fit with zero constraints on coefficients, such that lik is
	#	not significantly different
 	#    Inputs: 
	#	data: T x N  matrix of N  time series of length T
	#	p.mdl: selected order of VAR(p.mdl) model
	#	N.final: indices of final variables
 	#	alpha: p-value threshold for GLR test, for zero restrictions
	#	thresh: small epsilon value, like 10^(-12), used for help conditioning
	#    Outputs:
	#	obj: value of the Whittle likelihood
	#	phi.arry: array of N x N x p.mdl of the fitted VAR coefficients
	#	sigma: innovation covariance matrix of VAR
	#    Requires: package "fractal"
	#
	########################################################################		

data <- data[,N.final] 
T <- dim(data)[1]
N <- dim(data)[2]
mu <- colMeans(data)
data <- t(t(data) - mu)
data.acf <- acf(data,type="covariance",plot=FALSE,lag.max=T)$acf

## get basic unconstrained VAR fit, via YW
if(p.mdl==1) 
{
	data.toep <- data.acf[1,,]
	data.seg <- data.acf[2,,]
}
if(p.mdl==2) 
{
	data.toep <- rbind(cbind(data.acf[1,,],data.acf[2,,]),
		cbind(t(data.acf[2,,]),data.acf[1,,]))
	data.seg <- cbind(data.acf[2,,],data.acf[3,,])
}
if(p.mdl==3) 
{
	data.toep <- rbind(cbind(data.acf[1,,],data.acf[2,,],data.acf[3,,]),
			cbind(t(data.acf[2,,]),data.acf[1,,],data.acf[2,,]),
			cbind(t(data.acf[3,,]),t(data.acf[2,,]),data.acf[1,,]))
	data.seg <- cbind(data.acf[2,,],data.acf[3,,],data.acf[4,,])
}
if(p.mdl==4) 
{
	data.toep <- rbind(cbind(data.acf[1,,],data.acf[2,,],data.acf[3,,],data.acf[4,,]),
			cbind(t(data.acf[2,,]),data.acf[1,,],data.acf[2,,],data.acf[3,,]),
			cbind(t(data.acf[3,,]),t(data.acf[2,,]),data.acf[1,,],data.acf[2,,]),
			cbind(t(data.acf[4,,]),t(data.acf[3,,]),t(data.acf[2,,]),data.acf[1,,])) 
	data.seg <- cbind(data.acf[2,,],data.acf[3,,],data.acf[4,,],data.acf[5,,])
}
fit.yw <- ar.yw(data,aic=FALSE,order.max=p.mdl) 
phi.mat <- NULL
for(i in 1:p.mdl) 
{ 
	phi.mat <- cbind(phi.mat,fit.yw$ar[i,,])
}
sigma.hat <- data.acf[1,,] - phi.mat %*% data.toep %*% t(phi.mat)
gamma.acf <- acf(fit.yw$resid,na.action=na.pass,lag.max=40,plot=FALSE)
resids <- fit.yw$resid
obj <- T*log(det(sigma.hat))

# determine t statistics for coefficients
phi.arry <- aperm(fit.yw$ar,c(2,3,1))
asymp.var <- solve(data.toep) %x% sigma.hat
asymp.se <- sqrt(diag(asymp.var))/sqrt(T)
asymp.se <- array(matrix(asymp.se,nrow=N),c(N,N,p.mdl))
tstats <- NULL
for(i in 1:p.mdl) 
{ 
	tstats <- cbind(tstats,phi.arry[,,i]/asymp.se[,,i]) 
}
tstats <- matrix(tstats,ncol=1)
quants <- sort(abs(as.vector(tstats)))

# cycle through constrained YW models, until lik changes significantly
better <- TRUE
i <- 0
while(better) {
	i <- i+1
	quant <- quants[i]
	select.mat <- matrix(0,N,p.mdl*N)
	select.mat[abs(tstats) > quant] <- 1
	J.mat <- diag(p.mdl*N^2)
	ind <- seq(1,p.mdl*N^2)
	J.mat <- J.mat[,ind[as.vector(select.mat)==1]]
	r <- dim(J.mat)[2]
	#print(r)

	# Iterative method to solve
	psi <- matrix(phi.mat[select.mat==1],nrow=r,ncol=1)
	phi <- matrix(J.mat %*% psi,nrow=N,ncol=p.mdl*N)
	sigma <- data.acf[1,,] - phi %*% data.toep %*% t(phi)
	eps <- 1
#	thresh <- 10^(-12)
#	thresh <- 10^(-10)
	while(eps > thresh)
	{
		phi.old <- phi
		scale.mat <- t(J.mat) %*% (data.toep %x% solve(sigma)) %*% J.mat
#		scale.mat <- scale.mat + thresh*diag(r)
		psi <- solve(scale.mat) %*% t(J.mat) %*% (t(data.seg) %x% solve(sigma)) %*% 
			matrix(diag(N),nrow=N^2,ncol=1)
		phi <- matrix(J.mat %*% psi,nrow=N,ncol=p.mdl*N)
		sigma <- data.acf[1,,] - phi %*% t(data.seg) - data.seg %*% t(phi) + phi %*% data.toep %*% t(phi)
#		sigma <- sigma + thresh*diag(N)
		diff <- as.vector(matrix(phi - phi.old,nrow=p.mdl*N^2,ncol=1))
		eps <- sum(diff^2)
		#print(eps)
	}
	phi.second <- phi
	sigma.second <- sigma
	#print(phi.second)
 	#print(sigma.second)
	ldet.sigma <- sum(log(eigen(sigma.second)$values))
	print(c(T*ldet.sigma-obj,qchisq(1-alpha,df=(p.mdl*N^2-r))))
	better <- T*ldet.sigma-obj < qchisq(1-alpha,df=(p.mdl*N^2-r))
}

# get best restricted model
phi <- phi.old
phi.arry <- array(phi,c(N,N,p.mdl))
sigma <- data.acf[1,,]  - phi %*% t(data.seg) - data.seg %*% t(phi) + phi %*% data.toep %*% t(phi)
ldet.sigma <- sum(log(eigen(sigma)$values))
obj <- T*ldet.sigma

return(list(obj,phi.arry,sigma))
}



