specFactmvar <- function(xAcf)
{

	Nseries <- dim(xAcf)[1]
	q <- dim(xAcf)[2] - 1
	Zzero <- matrix(0,nrow=Nseries,ncol=Nseries)
	eps <- 1
	thresh <- 10^(-15)
	Linvt <- diag(Nseries)
	Dinv <- solve(xAcf[,1,])
	sqrtLam <- svd(xAcf[,1,])
	sqrtLam <- sqrtLam$u %*% diag(sqrt(sqrtLam$d)) %*% t(sqrtLam$u) 
	Gseq <- solve(t(sqrtLam))
	oldLam <- xAcf[,1,]
	gamSeq <- NULL
	count <- 1
	while(eps > thresh)
	{
		if(count <= q) nextgam <- xAcf[,count+1,] else nextgam <- Zzero
		gamSeq <- cbind(nextgam,gamSeq)
		Bseq <- gamSeq %*% Gseq
		Lam <- xAcf[,1,] - Bseq %*% t(Bseq)
		sqrtLam <- svd(Lam)
		sqrtLam <- sqrtLam$u %*% diag(sqrt(sqrtLam$d)) %*% t(sqrtLam$u) 
		Gseq <- rbind(cbind(Gseq,-1*Gseq%*%t(Bseq)%*%solve(t(sqrtLam))),
			cbind(t(rep(1,count) %x% Zzero),solve(t(sqrtLam))))
		count <- count+1
		if(count > q) eps <- sum((oldLam - Lam)^2)
		oldLam <- Lam
#		print(c(count,eps))
#		print(Lam)
#		print(max(eigen(Gseq)$values))
	}
	Thetas <- cbind(Bseq,sqrtLam) %*% (diag(count) %x% solve(sqrtLam))
	arrayThetas <- array(Thetas,dim=c(Nseries,Nseries,count))
	out <- list(arrayThetas[,,(count-q):count],Lam)
	return(out)
}
