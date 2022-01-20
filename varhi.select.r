varhi.select <- function(data,N.core,p.mdl,thresh)
{
  
  ########################################################################
  #	varhi.select by Tucker McElroy, March 19, 2021
  #
  #    Purpose: given a VAR(p) model, select variables
  #    Inputs: 
  #	data: T x N  matrix of N time series of length T,
  #		where  T corresponds to a subspan for forecast evaluation
  #	N.core: index subset of {1,2,...,N} corresponding to core series
  #	p.mdl: selected order of VAR(p.mdl) model
  # thresh: alpha value controlling FWER and type I error
  #    Outputs:
  # N.final: indices of final selected variables
  #
  ########################################################################		
  
  N <- dim(data)[2]
  wald.tests <- varhi.getaux(data,N.core,p.mdl)[[1]]
  ordering <- setdiff(seq(1,N),N.core)
  ordering <- rev(ordering[order(wald.tests)])
  N.build <- varhi.build(data,N.core,p.mdl,ordering,thresh)
  N.ult <- N.core
  if(length(N.build)>1)
  {
    pval <- varhi.refine(data[,N.build,drop=FALSE],N.core,p.mdl)[[2]]
    N.test <- setdiff(N.build,N.core)
    N.ult <- c(N.ult,N.test[pval < thresh])
  }
  N.final <- N.ult
 
  return(N.final)
}