getGUCD <- function(Sigma,Rank)
{

	#############################
	#	getGUCD
	#		by Tucker McElroy
	#	Gets Generalized Cholesky Decomposition of Sigma,
	#	 allowing for zero Schur complements,
	#	 in terms of upper triangular matrices.
	#      Rank is the presumed rank of the matrix, less than or equal
	#	  to the dimension of the matrix
	#	Output consists of the upper cholesky matrix,
	#	  and the diagonal matrix, of reduced dimension
	#
	#############################

	N <- dim(Sigma)[1]
	U.mat <- matrix(1,1,1)
	U.mat.inv <- U.mat
	C.mat <- Sigma[N,N]
	if(N > 1) {
	for(j in 2:N)
	{
		
		C.inv <- 1/C.mat
		C.inv[C.mat==0] <- 0
		new.sigma <- Sigma[(N+1-j),(N+2-j):N]
		if(j==2) { new.sigma <- as.matrix(new.sigma); U.mat <- as.matrix(U.mat) }
		new.u <- new.sigma %*% t(U.mat.inv)*C.inv
		new.u.tilde <- new.u %*% U.mat.inv
		U.mat <- cbind(rep(0,(j-1)),U.mat)
		U.mat <- rbind(c(1,new.u),U.mat)
		U.mat.inv <- cbind(rep(0,j-1),U.mat.inv)
		U.mat.inv <- rbind(c(1,-1*new.u.tilde),U.mat.inv)
		if(j==2) new.c <- Sigma[N-1,N-1] - new.u^2*C.mat
		if(j > 2) new.c <- Sigma[N+1-j,N+1-j] - new.u %*% diag(C.mat) %*% t(new.u)
		if(new.c <= 0) { new.c <- 0 }
		C.mat <- c(new.c,C.mat)
	} }
	
	rank.index <- rank(C.mat,ties.method="first")
	dims <- seq(1,N)[rank.index > (N-Rank)]

	U.mat <- matrix(U.mat[,dims],nrow=N,ncol=length(dims))
	C.mat <- C.mat[dims]

#	print(Lmat)
#	print(Dmat)
#	print(Lmat %*% diag(Dmat,nrow=length(dims)) %*% t(Lmat))

	return(list(U.mat,C.mat))
#	return(list(U.mat,C.mat,dims))
}



