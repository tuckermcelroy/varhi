VARMAauto <- function(phi,theta,sigma,maxlag)
{
	
	##################################################################################
	#
	#  VARMAauto
	#	Copyright (2015) Tucker McElroy
	#
	#	This program is free software; you can redistribute it and/or
	#	modify it under the terms of the GNU General Public License
	#	as published by the Free Software Foundation; either version 2
	#	of the License, or (at your option) any later version.
	#
	#	This program is distributed in the hope that it will be useful,
	#	but WITHOUT ANY WARRANTY; without even the implied warranty of
	#	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	#	GNU General Public License for more details.
	#
	#	You should have received a copy of the GNU General Public License
	#	along with this program; if not, write to the Free Software
	#	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
	#
	##############################################################################
	
	#######################################################################
	#	DOCUMENTATION:	
	# Function computes autocovariances of VARMA (p,q) from lag zero
	#	to maxlag, with inputs phi and theta.
	#	(1 - phi[1]z ... - phi[p]z^p) X_t = (1 + theta[1]z ...+ theta[q]z^q) WN
	#  output: autocovariance string of length maxlag
	#  for absent AR or MA portions, pass in NULL
	#  phi and theta should be arrays of m x m matrices
	#	sigma should be an m x m matrix
	#  e.g. phi <- array(cbind(phi1,phi2,...,phip),c(m,m,p))
	################################################################

polymulMat <- function(amat,bmat)
{
	p <- dim(amat)[3]
	q <- dim(bmat)[3]
	m <- dim(amat)[2]
	amatd <- amat[,,p:1]
	if(m==1) amatd <- array(amatd,c(m,m,p))
	if(q > 1) amatd <- array(c(matrix(0,m,m*(q-1)),amatd),c(m,m,p+q-1))
	bigmat <- NULL
	for(i in 1:(p+q-1)) 
	{
		nextmat <- matrix(amatd[,,1:(p+q-1)],m,m*(p+q-1))
		bigmat <- rbind(nextmat,bigmat)
		amatd <- amatd[,,-1]
		amatd <- array(c(amatd,matrix(0,m,m)),c(m,m,p+q-1))
	}
	bigmat <- bigmat[,1:(m*q)]
	out <- bigmat %*% t(matrix(bmat[,,q:1],m,m*q))
	out <- array(out,c(m,p+q-1,m))
	temp <- NULL
	for(i in 1:(p+q-1))
	{
		temp <- cbind(temp,out[,i,])
	}
	out <- array(temp,c(m,m,p+q-1))

	return(out)
}

Kcommut <- function(vect,m,n)
{
	return(matrix(t(matrix(vect,nrow=m,ncol=n)),ncol=1))
}

m <- dim(sigma)[2]
p <- 0
q <- 0
if (length(phi) > 0) p <- dim(phi)[3]
if (length(theta) > 0) q <- dim(theta)[3]
Kmat <- apply(diag(m^2),1,Kcommut,m,m)

if (q == 0) { gamMA <- array(sigma,c(m,m,1)) } else 
{
	temp <- polymulMat(array(cbind(diag(m),matrix(theta,m,m*q)),c(m,m,q+1)),
		array(sigma,c(m,m,1)))
	gamMA <- polymulMat(temp,array(cbind(diag(m),matrix(theta,m,m*q)),c(m,m,q+1)))
}
gamMA <- gamMA[,,(q+1):(2*q+1)]
if(m==1) gamMA <- array(gamMA,c(m,m,q+1))
gamMAvec <- matrix(gamMA,m^2*(q+1),1)

if (p > 0) 
{
	Amat <- matrix(0,nrow=m^2*(p+1),ncol=m^2*(2*p+1))
	Amat <- array(Amat,c(m^2,p+1,m^2,2*p+1))
	Arow <- diag(m^2)
	for(i in 1:p)
	{
		Arow <- cbind(-1*diag(m) %x% phi[,,i],Arow)
	}
	for(i in 1:(p+1))
	{
		Amat[,i,,i:(i+p)] <- Arow
	}
	newA <- array(matrix(Amat[,1:(p+1),,1:p],m^2*(p+1),m^2*(p)),c(m^2,p+1,m^2,p))
	for(i in 1:(p+1))
	{
		for(j in 1:p)
		{
 			newA[,i,,j] <- newA[,i,,j] %*% Kmat
		}
	}
	Amat <- cbind(matrix(Amat[,,,p+1],m^2*(p+1),m^2),
			matrix(Amat[,,,(p+2):(2*p+1)],m^2*(p+1),m^2*(p)) + 
			matrix(newA[,,,p:1],m^2*(p+1),m^2*(p)))

	Bmat <- matrix(0,nrow=m^2*(q+1),ncol=m^2*(p+q+1))
	Bmat <- array(Bmat,c(m^2,q+1,m^2,p+q+1))
	Brow <- diag(m^2)
	for(i in 1:p)
	{
		Brow <- cbind(Brow,-1*phi[,,i] %x% diag(m))
	}
	for(i in 1:(q+1))
	{
		Bmat[,i,,i:(i+p)] <- Brow
	}
	Bmat <- Bmat[,,,1:(q+1)]
	Bmat <- matrix(Bmat,m^2*(q+1),m^2*(q+1))
	Binv <- solve(Bmat)

	gamMix <- Binv %*% gamMAvec
	if (p <= q) gamMixTemp <- gamMix[1:((p+1)*m^2)] else 
		gamMixTemp <- c(gamMix,rep(0,(p-q)*m^2))
	gamARMA <- solve(Amat) %*% gamMixTemp 
	gamMix <- array(matrix(gamMix,m,m*(q+1)),c(m,m,q+1))
	gamARMA <- array(matrix(gamARMA,m,m*(p+1)),c(m,m,p+1))
} else 
{
	gamARMA <- array(gamMA[,,1],c(m,m,1))
	if (q == 0) { gamMix <- array(sigma,c(m,m,1)) } else 	
		{ 
			gamMix <- gamMA[,,1:(q+1)]
			if(m==1) gamMix <- array(gamMix,c(1,1,q+1)) 
		}
}

if (maxlag <= p) 
{
	gamARMA <- gamARMA[,,1:(maxlag+1)] 
	if(m==1) gamARMA <- array(gamARMA,c(1,1,maxlag+1)) 
} else
{
	if (maxlag > q) gamMix <- array(cbind(matrix(gamMix,m,m*(q+1)),
		matrix(0,m,m*(maxlag-q))),c(m,m,(maxlag+1)))
	for(k in 1:(maxlag-p))
	{
		len <- dim(gamARMA)[3]
		acf <- gamMix[,,p+1+k]
		if (p > 0) 
		{
			temp <- NULL
			for(i in 1:p)
			{
				temp <- rbind(temp,gamARMA[,,len-i+1])
			}
			acf <- acf + matrix(phi,m,m*p) %*% temp
		} 
		gamARMA <- array(cbind(matrix(gamARMA,m,m*len),acf),c(m,m,len+1))
	}
}

return(gamARMA)
}
