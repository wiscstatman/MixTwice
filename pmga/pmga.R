
# Some code to implement part of the Pool Monotone groups
# algorithm.  Started 12/14/97

pmga <- function(ss,ff)
	{
	# the pool monotone groups algorithm; from MAN 12/14/97

	# ss,ff are vectors of equal lengths, non-negative, ss+ff>0
	# representing successes and failures at a set of sample points.
	# aiming to estimate MLE of success probs subject to these probs being non-decreasing

	nn <- length(ss)  # I assume nn > 1
	prop <- ss/(ss+ff)
	ind <- numeric(nn)
	ind[1] <- 1
	for( j in 2:nn ) # identify non-increasing groups
	 {
	  ind[j] <- ifelse( prop[j] <= prop[j-1], ind[j-1], ind[j-1]+1 )
	 }
	tracer <- ind
	mm <- ind[nn]
	if( mm==nn )
	 {
	  ss2 <- ss
	  ff2 <- ff
	 }
	else # pool groups
	 {
	  ss2 <- numeric(mm)
	  ff2 <- numeric(mm)
	  for( j in 1:mm )
	   {
	    ss2[j] <- sum( ss[ind==j] )
	    ff2[j] <- sum( ff[ind==j] )
	    tracer[ind==j] <- j
	   }
	 }
	return( list( ss=ss2, ff=ff2, tracer=tracer ) )
	}


# try it

set.seed(453)

x <- rbinom(10, prob=seq(.1,.9, length=10), size=20)
f <- pmga(x,20-x)

fit <- f$ss[f$tracer]/( f$ss[f$tracer] + f$ff[f$tracer] )
