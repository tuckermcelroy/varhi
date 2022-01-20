spec.ridge <- function(gamma,eps,grid)
{

  #############################################
  #  spec.ridge by Tucker McElroy, 1/24/2020
  #
  #  Purpose: compute a ridged spectral density estimator 
  #  Inputs:
  #   gamma: N x N x (q+1) array of autocovariances
  #   eps: a small non-negative number, for numerical stability
  #   grid: mesh size for frequencies
  #  Outputs:
  #   f.ridge: N x N x (grid+1) spectral density estimate
  #   acf.ridge: N x N x (q+1) array of ridged autocovariances
  #
  ############################################

  N <- dim(gamma)[1]
  q <- dim(gamma)[3]-1
  gam0.min <- min(eigen(gamma[,,1])$value)
  lambda <- pi*seq(0,grid)/grid

  spec.f <- t(exp(-1i*lambda*0)) %x% gamma[,,1]
  for(k in 1:q)
  {
    spec.f <- spec.f + t(exp(-1i*lambda*k)) %x% gamma[,,k+1] + t(exp(1i*lambda*k)) %x% t(gamma[,,k+1])
  }
  f.array <- array(spec.f,c(N,N,grid+1))
  f.min <- Inf
  for(j in 0:grid)
  {
     f.min <- min(f.min,min(Re(eigen(f.array[,,j+1])$value)))
  }
  if(f.min > 0) { alpha <- 1 } else { alpha <- gam0.min/(gam0.min - f.min) - eps }
  f.ridge <- alpha*spec.f + (1-alpha)*t(exp(-1i*lambda*0)) %x% gamma[,,1]
  f.ridge <- array(f.ridge,c(N,N,grid+1)) 
  
#  spec.g <- t(0*exp(-1i*lambda*0)) %x% gamma[,,1]
#  for(k in 1:q)
#  {
#    spec.g <- spec.g + t(exp(-1i*lambda*k)) %x% gamma[,,k+1] + t(exp(1i*lambda*k)) %x% t(gamma[,,k+1])
#  }
#  spec.f <- spec.g
#  spec.g <- array(spec.g,c(N,N,grid+1))
#  g.min <- Inf
#  for(j in 0:grid)
#  {
#     g.min <- min(g.min,min(eigen(spec.g[,,j+1])$value))
#  }
#  if(g.min > 0) { alpha <- 1 } else { alpha <- -1*gam0.min/g.min - eps }
#  spec.f <- alpha*spec.f + t(exp(-1i*lambda*0)) %x% gamma[,,1]

  acf.ridge <- array(0,c(N,N,q+1))
  acf.ridge[,,1] <- gamma[,,1]
  for(k in 1:q) { acf.ridge[,,k+1] <- alpha*gamma[,,k+1] }
#  print(alpha)
    
  return(list(f.ridge,acf.ridge))
}
