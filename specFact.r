specFact <- function(poly)
{
	# input a 2-sided polynomial (symmetric coefficients)
	#  output is a polynomial theta, which should have the property that
	#   polymult(theta,rev(theta)) = poly
	p0 <- length(poly)-1
	prod <- poly[1]
	poly <- poly/poly[1]
	roots <- polyroot(poly)
	theta <- 1
	for(i in 1:p0)
		if (Mod(roots[i]) > 1)
		{
#			theta <- polymult(theta,c(1,-1/roots[i]))
			theta <- polymult(theta,c(-1*roots[i],1))
			prod <- prod/(-1*roots[i])
#			print(c(i,theta))
#			print(Mod(polyroot(theta)))
		}
	theta <- theta*sqrt(prod)
	return(theta)
}
			