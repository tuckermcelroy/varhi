varhi.genvar <- function(phi,N.core,N.aux,N.anc,p.mdl)
{
  
  ########################################################################
  #	varhi.genvar by Tucker McElroy, February 17, 2021
  #
  #    Purpose: given variables and a VAR model order,  
  #     return a stable VAR matrix, where the upper right block
  #     corresponding to ancillary variables is zero
  #    Inputs: 
  # phi: p.mdl*N^2 dimensional vector of initial values, where
  #   N is composed of N.core + N.aux + N.anc
  #	N.core: index subset of {1,2,...,N} corresponding to core series
  #	N.aux: index subset of {1,2,...,N} corresponding to auxiliary series
  #	N.anc: index subset of {1,2,...,N} corresponding to ancillary series
  #	p.mdl: selected order of VAR(p.mdl) model
  #    Outputs:
  #	phi.stable: array N x N x p.mdl of stable VAR coefficient matrices
  #
  #  Requires: VARMAauto, specFactmvar
  ########################################################################		
  
  var.stabilize <- function(phi.prime)
  {
    N <- dim(phi.prime)[1]
    p <- dim(phi.prime)[2]/N
    phi.array <- array(-1*t(phi.prime),c(N,N,p))
    gamma.inverse <- VARMAauto(NULL,phi.array,diag(N),p)
    gamma.inverse.alt <- array(0,c(N,p+1,N))
    for(i in 1:(p+1)) { gamma.inverse.alt[,i,] <- gamma.inverse[,,i] }
    out <- specFactmvar(gamma.inverse.alt)
    if(p > 0) phi.inverse.array <- out[[1]][,,p:1]
    if(p==1) phi.inverse.array <- array(out[[1]][,,1],c(N,N,1))
    phi.stable.array <- array(0,c(N,N,p))
    for(i in 1:p) { phi.stable.array[,,i] <- -1*t(phi.inverse.array[,,i]) }
    phi.stable <- matrix(phi.stable.array,c(N,p*N))
    return(phi.stable.array)
  }

  N.null <- union(N.core,N.aux)
  N.all <- union(N.null,N.anc)
  N <- length(N.all)
  phi.prime <- matrix(phi,N,p.mdl*N)
  sigma.prime <- diag(N)
  phi.prime.array <- array(phi.prime,c(N,N,p.mdl))
  phi.stable <- phi.prime.array
  phi.null.array <- phi.prime.array[N.null,N.null,,drop=FALSE] 
  phi.null.stable <- var.stabilize(matrix(phi.null.array,c(length(N.null),p.mdl*length(N.null))))
  for(i in 1:p.mdl)
  {
     phi.stable[N.null,N.null,i] <- phi.null.stable[,,i]
   }
  if(length(N.anc)>0) 
  { 
    phi.anc.array <- phi.prime.array[N.anc,N.anc,,drop=FALSE] 
    phi.anc.stable <- var.stabilize(matrix(phi.anc.array,c(length(N.anc),p.mdl*length(N.anc))))
    for(i in 1:p.mdl)
    {
      phi.stable[N.anc,N.anc,i] <- phi.anc.stable[,,i]
      phi.stable[N.null,N.anc,i] <- 0
    }
  }
  return(phi.stable)  
} 
 
  