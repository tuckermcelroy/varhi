varhi.lasso <- function(data.sub,p.mdl,lambdas,ywinit=FALSE)
{
  
  ########################################################################
  #	varhi.lasso by Tucker McElroy, January 27, 2021
  #
  #    Purpose: fit VAR lasso model to data.
  #	We fit the core+auxiliary model, with group lasso penalty for
  # each row and column of all coefficient matrices
  #
  #    Inputs: 
  #
  #	data.sub: T x N  matrix of N time series of length T,
  #		where  T corresponds to a subspan for forecast evaluation
  #	p.mdl: selected order of VAR(p.mdl) model
  # lambdas: LASSO penalties, one for each variable (N of them).
  #   To retain a core variable, set corresponding lambda to zero.
  # ywinit: set to TRUE if the optimization should be initialized with
  #   the Yule-Walker estimate. Defaults to FALSE, where initialization is 0.
  #
  #    Outputs: 
  #
  # list(var.set,phi.array,fev), where var.set is subset of 1,...,N
  #   of selected variables.  
  #   phi.array and fev are corresponding VAR coefficients and Sigma
  #
  #
  ########################################################################		
  
  fev_mat <- function(phi.array,gamma)
  {
    p.mdl <- dim(phi.array)[3]
    N <- dim(phi.array)[2]
    phi.mat <- matrix(phi.array,c(N,p.mdl*N))
    fev <- gamma[,,1] 
    for(j in 1:p.mdl)
    { 
      fev <- fev - phi.array[,,j] %*% t(gamma[,,j+1]) - gamma[,,j+1] %*% t(phi.array[,,j])
      for(k in 1:p.mdl)
      {
        lag <- k-j
        if(lag >= 0) { gamma_mat <- gamma[,,lag+1] } else { gamma_mat <- t(gamma[,,-lag+1]) }
        fev <- fev + phi.array[,,j] %*% gamma_mat %*% t(phi.array[,,k])
      }
    }
    return(fev)
  }
  
  lasso_crit <- function(phi,gamma,lambdas)
  {
    N <- dim(gamma)[1]
    p.mdl <- length(phi)/N^2
    phi.array <- array(0,c(N,N,p.mdl))
    con <- rep(0,N)
    for(k in 1:p.mdl)
    {
      phi.array[,,k] <- matrix(phi[((k-1)*N^2 + 1):(k*N^2)],c(N,N))
      for(i in 1:N)
      {
        con[i] <- con[i] + sum(abs(phi.array[,i,k])) 
        #+ sum(abs(phi.array[i,,k]))
      }
    }
    fev <- fev_mat(phi.array,gamma)
    lik <- log(det(fev)) + sum(lambdas * con)
    return(lik)
  }

  thresh <- 10^{-12}
  x <- t(data.sub)
  N <- dim(x)[1]
  T <- dim(x)[2]
  x.acf <- acf(t(x),plot=FALSE,lag.max = p.mdl,type="covariance")$acf
  gamma <- aperm(x.acf,c(2,3,1))
  gam0.min <- min(eigen(gamma[,,1])$value)
  eps <- 10^(-6)
  gamma[,,1] <- gamma[,,1] + (gam0.min + eps)* diag(N)

  phi.init <- rep(0,p.mdl*N^2)
  if(ywinit) phi.init <- matrix(matrix(aperm(ar.yw(t(x),order=p.mdl,aic=FALSE)$ar,c(2,3,1)),nrow=N),ncol=1)
  lasso.fit <- optim(phi.init,lasso_crit,gamma=gamma,lambdas=lambdas,method="BFGS")  
#  print(lasso.fit)
  
  phi.lasso <- lasso.fit$par
  phi.array <- array(0,c(N,N,p.mdl))
  con <- rep(0,N)
  for(k in 1:p.mdl)
  {
    phi.array[,,k] <- matrix(phi.lasso[((k-1)*N^2 + 1):(k*N^2)],c(N,N))
    for(i in 1:N)
    {
      con[i] <- con[i] + sum(abs(phi.array[,i,k])) 
      #+ sum(abs(phi.array[i,,k]))
    }
  }
  fev <- fev_mat(phi.array,gamma)
  
  var.set <- seq(1,N)[con > thresh]
  phi.array <- phi.array[var.set,var.set,,drop=FALSE]
  fev <- fev[var.set,var.set,drop=FALSE]
  
  return(list(var.set,phi.array,fev))
}