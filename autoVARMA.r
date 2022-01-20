autoVARMA <- function(phi,theta,phiseas,thetaseas,sigma,season,grid=1000,maxlag)
{
  
  ##########################################################################
  #
  #	autoVARMA
  # 	    Copyright (C) 2020  Tucker McElroy
  #
  #    This program is free software: you can redistribute it and/or modify
  #    it under the terms of the GNU General Public License as published by
  #    the Free Software Foundation, either version 3 of the License, or
  #    (at your option) any later version.
  #
  #    This program is distributed in the hope that it will be useful,
  #    but WITHOUT ANY WARRANTY; without even the implied warranty of
  #    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  #    GNU General Public License for more details.
  #
  #    You should have received a copy of the GNU General Public License
  #    along with this program.  If not, see <https://www.gnu.org/licenses/>.
  #
  ############################################################################
  
  ################# Documentation #####################################
  #
  #	Purpose: computes autocovariances of SVARMA usng frequency domain
  #	Background: function computes autocovariances of SVARMA (p,q,ps,qs) from lag zero
  #		to maxlag, with array inputs phi and theta.  VARMA equation:
  #	(1 - phi[1]B ... - phi[p]B^p) (1 - Phi[1]B^s ... - Phi[ps]B^{s*ps}) X_t =
  #    (1 + theta[1]B ... + theta[q]B^q) (1 + Theta[1]B^s ... + Theta[qs]B^{s*qs}) WN_t.
  #	Note: for absent VAR or VMA portions, pass in NULL
  #	Inputs:
  #		phi: array of dimension m x m x p of VAR coefficients, e.g.,
  #			phi <- array(cbind(phi1,phi2,...,phip),c(m,m,p))
  #		theta: array of dimension m x m x q of VMA coefficients, e.g.,
  #			theta <- array(cbind(theta1,theta2,...,thetaq),c(m,m,q))
  #		phiseas: array of dimension m x m x ps of SVAR coefficients, e.g.,
  #			phi <- array(cbind(phi1,phi2,...,phips),c(m,m,ps))
  #		thetaseas: array of dimension m x m x qs of SVMA coefficients, e.g.,
  #			theta <- array(cbind(theta1,theta2,...,thetaqs),c(m,m,qs))
  #		sigma: m x m covariance matrix of white noise
  #   grid: Riemann mesh size
  #	Outputs:
  #		autocovariances at lags 0 through maxlag, as array of dimension m x m x (maxlag+1)
  #
  ####################################################################
  
  polymult <- function(a,b)
  {
    bb <- c(b,rep(0,length(a)-1))
    B <- toeplitz(bb)
    B[lower.tri(B)] <- 0
    aa <- rev(c(a,rep(0,length(b)-1)))
    prod <- B %*% matrix(aa,length(aa),1)
    return(rev(prod[,1]))
  }
  
  polymulMat <- function(amat,bmat)
  {
    p <- dim(amat)[3]-1
    q <- dim(bmat)[3]-1
    N <- dim(amat)[2]
    
    r <- p+q
    bmat.pad <- array(0,c(N,N,r+1))
    for(i in 1:(q+1)) { bmat.pad[,,i] <- bmat[,,i] }
    cmat <- array(0,c(N,N,r+1))
    cmat[,,1] <- amat[,,1] %*% bmat.pad[,,1]
    for(j in 2:(r+1))
    {
      cmat[,,j] <- amat[,,1] %*% bmat.pad[,,j]
      for(k in 1:min(p,j-1))
      { cmat[,,j] <- cmat[,,j] + amat[,,k+1] %*% bmat.pad[,,j-k] }
    }
    
    return(cmat)
  }
  
  ar.adjoint <- function(poly.array)
  {
    p <- dim(poly.array)[3]-1
    N <- dim(poly.array)[1]
    poly.0 <- poly.array[,,1]
    ar.poly <- det(poly.0)
    r <- p*(N-1)+1
    adj.array <- array(0,c(N,N,r))
    adj.array[,,1] <- ar.poly*solve(poly.0)
    
    if(p > 0) {
      poly.mat <- matrix(poly.array[,,2:(p+1)],c(N,p*N))
      poly.coefs <- -1*solve(poly.0) %*% poly.mat
      
      if(p==1) { companion.mat <- poly.coefs } else {
        companion.mat <- diag(p*N)[1:((p-1)*N),]
        companion.mat <- rbind(poly.coefs,companion.mat) }
      poly.evals <- eigen(companion.mat)$values
      for(j in 1:(p*N))
      { ar.poly <- polymult(ar.poly,c(1,-poly.evals[j])) }
      ar.poly <- Re(ar.poly)
      
      for(j in 2:r)
      {
        adj.array[,,j] <- ar.poly[j]*solve(poly.0)
        for(k in 1:min(p,j-1))
        {
          adj.array[,,j] <- adj.array[,,j] - solve(poly.0) %*%
            poly.array[,,k+1] %*% adj.array[,,j-k]
        }
      } }
    
    return(list(adj.array,ar.poly))
  }
  
  N <- dim(sigma)[1]
  p <- 0
  q <- 0
  ps <- 0
  qs <- 0
  if (length(phi) > 0) p <- dim(phi)[3]
  if (length(theta) > 0) q <- dim(theta)[3]
  if (length(phiseas) > 0) ps <- dim(phiseas)[3]
  if (length(thetaseas) > 0) qs <- dim(thetaseas)[3]
  phi.long <- array(diag(N),c(N,N,1))
  if(p > 0) { phi.long <- array(cbind(diag(N),-1*matrix(phi,c(N,N*p))),c(N,N,p+1)) }
  out <- ar.adjoint(phi.long)
  phi.adjoint <- out[[1]]
  phi.det <- out[[2]]
  theta.long <- array(diag(N),c(N,N,1))
  if(q > 0) { theta.long <- array(cbind(diag(N),matrix(theta,c(N,N*q))),c(N,N,q+1)) }
  phiseas.long <- array(diag(N),c(N,N,1))
  if(ps > 0) { phiseas.long <- array(cbind(diag(N),-1*matrix(phiseas,c(N,N*ps))),c(N,N,ps+1)) }
  out <- ar.adjoint(phiseas.long)
  phiseas.adjoint <- out[[1]]
  phiseas.det <- out[[2]]
  thetaseas.long <- array(diag(N),c(N,N,1))
  if(qs > 0) { thetaseas.long <- array(cbind(diag(N),matrix(thetaseas,c(N,N*qs))),c(N,N,qs+1)) }
  freqs <- seq(-1*grid,grid)*pi/grid
  len <- length(freqs)
  lambdas <- exp(-1i*freqs)
  
  phi.adjointz <- array(0,c(N,N*len))
  for(k in 1:dim(phi.adjoint)[3])
  {
    phi.adjointz <- phi.adjointz + t(lambdas^(k-1)) %x% phi.adjoint[,,k]
  }
  phi.adjointz <- array(phi.adjointz,c(N,N,len))
  phi.detz <- rep(0,len)
  for(k in 1:length(phi.det))
  {
    phi.detz <- phi.detz + lambdas^(k-1) * phi.det[k]
  }
  theta.z <- array(0,c(N,N*len))
  for(k in 1:dim(theta.long)[3])
  {
    theta.z <- theta.z + t(lambdas^(k-1)) %x% theta.long[,,k]
  }
  theta.z <- array(theta.z,c(N,N,len))
  phiseas.adjointz <- array(0,c(N,N*len))
  for(k in 1:dim(phiseas.adjoint)[3])
  {
    phiseas.adjointz <- phiseas.adjointz + t(lambdas^(season*(k-1))) %x% phiseas.adjoint[,,k]
  }
  phiseas.adjointz <- array(phiseas.adjointz,c(N,N,len))
  phiseas.detz <- rep(0,len)
  for(k in 1:length(phiseas.det))
  {
    phiseas.detz <- phiseas.detz + lambdas^(season*(k-1)) * phiseas.det[k]
  }
  thetaseas.z <- array(0,c(N,N*len))
  for(k in 1:dim(thetaseas.long)[3])
  {
    thetaseas.z <- thetaseas.z + t(lambdas^(season*(k-1))) %x% thetaseas.long[,,k]
  }
  thetaseas.z <- array(thetaseas.z,c(N,N,len))
  
  gamma <- array(0,c(N,N,maxlag+1))
  for(h in 0:maxlag)
  {
    val <- do.call(cbind,lapply(seq(1,len),function(i)
      phiseas.adjointz[,,i] %*% phi.adjointz[,,i] %*%
        theta.z[,,i] %*% thetaseas.z[,,i] %*% sigma %*%
        t(Conj(thetaseas.z[,,i])) %*% t(Conj(theta.z[,,i])) %*%
        t(Conj(phi.adjointz[,,i])) %*% t(Conj(phiseas.adjointz[,,i])) *
        Mod(phi.detz[i])^{-2} * Mod(phiseas.detz[i])^{-2} * lambdas[i]^{-h}))
    val <- len^{-1}*val %*% (rep(1,len) %x% diag(N))
    gamma[,,h+1] <- Re(val)
  }
  
  return(gamma)
  
}