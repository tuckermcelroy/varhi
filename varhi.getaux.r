varhi.getaux <- function(data.sub,N.core,p.mdl)
{

	########################################################################
	#	varhi.getaux by Tucker McElroy, August 15, 2017
	#
	#    Purpose: given core variables and a VAR model order, compute the
	#	impact of each ancillary variable, according to wald tests.
	#	We test the core+auxiliary model versus the core model, for each
	#	auxiliary variable
 	#    Inputs: 
	#	data.sub: T x N  matrix of N time series of length T,
	#		where  T corresponds to a subspan for forecast evaluation
	#	N.core: index subset of {1,2,...,N} corresponding to core series
	#	p.mdl: selected order of VAR(p.mdl) model
	#    Outputs:
	#	wald.tests: vector of length N giving all wald statistics	
	#
	########################################################################		

x <- t(data.sub)
N <- dim(x)[1]
T <- dim(x)[2]

N.anc <- setdiff(seq(1,N),N.core)
wald.tests <- NULL
pvals <- NULL

for(N.aux in N.anc) {

# fit models consisting of core plus an auxiliary variable

N.all <- union(N.core,N.aux)
data <- as.matrix(t(x)[,N.all])
mu <- colMeans(data)
data <- t(t(data) - mu)
data.acf <- acf(data,type="covariance",plot=FALSE)$acf
if(p.mdl==1) { 
	fit.order1 <- ar.yw(data,aic=FALSE,order.max=1) 
	if(dim(data)[2]>1) {
		phi.hat <- fit.order1$ar[1,,]
		sigma.hat <- data.acf[1,,] - phi.hat %*% data.acf[1,,] %*% t(phi.hat)
	} else
	{
		phi.hat <- fit.order1$ar
		sigma.hat <- (1 - phi.hat^2) * data.acf[1,,]
	}
#	acf(fit.order1$resid,na.action=na.pass)
}
if(p.mdl==2) {
	fit.order2 <- ar.yw(data,aic=FALSE,order.max=2)
	if(dim(data)[2]>1) {
		phi1.hat <- fit.order2$ar[1,,]
		phi2.hat <- fit.order2$ar[2,,]
		sigma.hat <- data.acf[1,,] - cbind(phi1.hat,phi2.hat) %*% 
			rbind(cbind(data.acf[1,,],data.acf[2,,]),cbind(t(data.acf[2,,]),data.acf[1,,])) %*% 
			t(cbind(phi1.hat,phi2.hat))
	} else
	{
		phi1.hat <- fit.order2$ar[1]
		phi2.hat <- fit.order2$ar[2]
		sigma.hat <- data.acf[1] - matrix(c(phi1.hat,phi2.hat),1,2) %*% 
			matrix(c(data.acf[1,,],data.acf[2,,],data.acf[2,,],data.acf[1,,]),2,2) %*% 
 				t(matrix(c(phi1.hat,phi2.hat),1,2))
	}
#	acf(fit.order2$resid,na.action=na.pass)
}
if(p.mdl==3) {
	fit.order3 <- ar.yw(data,aic=FALSE,order.max=3)
	if(dim(data)[2]>1) {
		phi1.hat <- fit.order3$ar[1,,]
		phi2.hat <- fit.order3$ar[2,,]
		phi3.hat <- fit.order3$ar[3,,]
		sigma.hat <- data.acf[1,,] - cbind(phi1.hat,phi2.hat,phi3.hat) %*% 
			rbind(cbind(data.acf[1,,],data.acf[2,,],data.acf[3,,]),
				cbind(t(data.acf[2,,]),data.acf[1,,],data.acf[2,,]),
				cbind(t(data.acf[3,,]),t(data.acf[2,,]),data.acf[1,,])) %*% 
			t(cbind(phi1.hat,phi2.hat,phi3.hat))
	} else
	{
		phi1.hat <- fit.order3$ar[1]
		phi2.hat <- fit.order3$ar[2]
		phi3.hat <- fit.order3$ar[3]
		sigma.hat <- data.acf[1] - matrix(c(phi1.hat,phi2.hat,phi3.hat),1,3) %*% 
			matrix(c(data.acf[1,,],data.acf[2,,],data.acf[3,,],data.acf[2,,],data.acf[1,,],
				data.acf[2,,],data.acf[3,,],data.acf[2,,],data.acf[1,,]),3,3) %*% 
 				t(matrix(c(phi1.hat,phi2.hat,phi3.hat),1,3))
	}
#	acf(fit.order3$resid,na.action=na.pass)
}
if(p.mdl==4) {
	fit.order4 <- ar.yw(data,aic=FALSE,order.max=4)
	if(dim(data)[2]>1) {
		phi1.hat <- fit.order4$ar[1,,]
		phi2.hat <- fit.order4$ar[2,,]
		phi3.hat <- fit.order4$ar[3,,]
		phi4.hat <- fit.order4$ar[4,,]
		sigma.hat <- data.acf[1,,] - cbind(phi1.hat,phi2.hat,phi3.hat,phi4.hat) %*% 
			rbind(cbind(data.acf[1,,],data.acf[2,,],data.acf[3,,],data.acf[4,,]),
				cbind(t(data.acf[2,,]),data.acf[1,,],data.acf[2,,],data.acf[3,,]),
				cbind(t(data.acf[3,,]),t(data.acf[2,,]),data.acf[1,,],data.acf[2,,]),
				cbind(t(data.acf[4,,]),t(data.acf[3,,]),t(data.acf[2,,]),data.acf[1,,])) %*% 
			t(cbind(phi1.hat,phi2.hat,phi3.hat,phi4.hat))
	} else
	{
		phi1.hat <- fit.order4$ar[1]
		phi2.hat <- fit.order4$ar[2]
		phi3.hat <- fit.order4$ar[3]
		phi4.hat <- fit.order4$ar[4]
		sigma.hat <- data.acf[1] - matrix(c(phi1.hat,phi2.hat,phi3.hat,phi4.hat),1,4) %*% 
			matrix(c(data.acf[1,,],data.acf[2,,],data.acf[3,,],data.acf[4,,],
				data.acf[2,,],data.acf[1,,],data.acf[2,,],data.acf[3,,],
				data.acf[3,,],data.acf[2,,],data.acf[1,,],data.acf[2,,],			
				data.acf[4,,],data.acf[3,,],data.acf[2,,],data.acf[1,,]),4,4) %*% 
 				t(matrix(c(phi1.hat,phi2.hat,phi3.hat,phi4.hat),1,4))
	}
#	acf(fit.order4$resid,na.action=na.pass)
}

## do Wald tests

E1 <- matrix(diag(length(N.all))[1:length(N.core),],ncol=length(N.all))
E2 <- matrix(diag(length(N.all))[(1+length(N.core)):length(N.all),],ncol=length(N.all))
E2 <- t(E2)
if(p.mdl==1) {
phi.core.aux <- matrix(E1 %*% phi.hat %*% E2,ncol=1)
var.core.aux <- t(E2) %*% solve(data.acf[1,,]) %*% E2 
var.core.aux <- var.core.aux %x% (E1 %*% sigma.hat %*% t(E1))
}
if(p.mdl==2) {
phi.core.aux <- matrix(cbind(E1 %*% phi1.hat %*% E2, E1 %*% phi2.hat %*% E2),ncol=1)
var.core.aux <- solve(rbind(cbind(data.acf[1,,],data.acf[2,,]),cbind(t(data.acf[2,,]),data.acf[1,,])))
var.core.aux <- (diag(p.mdl) %x% t(E2)) %*% var.core.aux %*% (diag(p.mdl) %x% E2) 
var.core.aux <- var.core.aux %x% (E1 %*% sigma.hat %*% t(E1))
}
if(p.mdl==3) {
phi.core.aux <- matrix(cbind(E1 %*% phi1.hat %*% E2, E1 %*% phi2.hat %*% E2, 
			E1 %*% phi3.hat %*% E2),ncol=1)
var.core.aux <- solve(rbind(cbind(data.acf[1,,],data.acf[2,,],data.acf[3,,]),
				cbind(t(data.acf[2,,]),data.acf[1,,],data.acf[2,,]),
				cbind(t(data.acf[3,,]),t(data.acf[2,,]),data.acf[1,,])))
var.core.aux <- (diag(p.mdl) %x% t(E2)) %*% var.core.aux %*% (diag(p.mdl) %x% E2) 
var.core.aux <- var.core.aux %x% (E1 %*% sigma.hat %*% t(E1))
}
if(p.mdl==4) {
phi.core.aux <- matrix(cbind(E1 %*% phi1.hat %*% E2, E1 %*% phi2.hat %*% E2, 
			E1 %*% phi3.hat %*% E2, E1 %*% phi4.hat %*% E2),ncol=1)
var.core.aux <- solve(rbind(cbind(data.acf[1,,],data.acf[2,,],data.acf[3,,],data.acf[4,,]),
				cbind(t(data.acf[2,,]),data.acf[1,,],data.acf[2,,],data.acf[3,,]),
				cbind(t(data.acf[3,,]),t(data.acf[2,,]),data.acf[1,,],data.acf[2,,]),
				cbind(t(data.acf[4,,]),t(data.acf[3,,]),t(data.acf[2,,]),data.acf[1,,])))
var.core.aux <- (diag(p.mdl) %x% t(E2)) %*% var.core.aux %*% (diag(p.mdl) %x% E2) 
var.core.aux <- var.core.aux %x% (E1 %*% sigma.hat %*% t(E1))
}

wald.test <- T* t(phi.core.aux) %*% solve(var.core.aux) %*% phi.core.aux
wald.tests <- c(wald.tests,wald.test)
pvals <- c(pvals,1-pchisq(wald.test,df=length(phi.core.aux)))

#print(c(wald.test,1-pchisq(wald.test,df=length(phi.core.aux))))

}


return(list(wald.tests,pvals))
}

