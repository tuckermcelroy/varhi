VAR2Wold <- function(phi,scale,maxlag)
{

	#################################
	#
	#  VAR2Wold by Tucker McElroy
	#
	#	input phi is array m x m x p
	#	output is array m x m x maxlag+1
	#	 of Wold coefficients, each right multiplied by scale matrix
	#
	################################

p <- 0
if (length(phi) > 0) p <- dim(phi)[3]
m <- dim(phi)[2]

wold <- array(0,c(m,m,maxlag+1))
wold[,,1] <- scale
for(i in 2:(maxlag+1))
{
	k <- min(i-1,p)
	wold.next <- 0*diag(m)
	for(l in 1:k)
	{
		wold.next <- wold.next + phi[,,l] %*% wold[,,i-l]
	}
	wold[,,i] <- wold.next 
}

return(wold)
}

