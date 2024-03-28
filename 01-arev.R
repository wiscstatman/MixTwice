
rm(list=ls())

load("peptideRAStats.RData")
library(fdrtool)

x <- array1_weak

tt <- x[,1]/sqrt(x[,2])

ff <- approxfun( density(tt) )


ff.val <- ff(tt)
f0.val <- dt(tt, df=22)

## from the p-value perspective

pv <- 2*pt( abs(tt), df=22, lower.tail=FALSE )
g <- fdrtool(pv, statistic="pvalue" )  ## 

## try to use that pi0 in a density estimate on tt scale 
## [note that locfdr gives pi0 = 1, and doesn't allow me to fix]

pi0 <- g$param[3]
fA <- ( ff.val - pi0*f0.val )/(1-pi0)  ## too few (i.e. neg density) at small negative tt

lfdr.uv <- pi0*f0.val/ff.val
lfdr.u <- g$lfdr

## 
fA.fun <- function(x, pi0=0.807)
 {
   tmp <- ( ff(x) - pi0*dt(x,df=22) )/(1-pi0)
   res <- ifelse( tmp >=0, tmp, 0 )
   res
 }


### there is something arbitrary about locfdr's choice of z-score scale
## for the RA data, locfdr can't get nulltype=0 to work
## but what if we worked with uniform scale

uu <- pt( tt, df=22 )

## hist(uu, prob=TRUE )  ## the fit will depend a bit on the choice of bins and mass at 0

nbin <- 201
br <- c( 0, seq( (min(uu)/2),  (max(uu)+1)/2, length=nbin-1 ) )
h <- hist(uu,  plot=FALSE, breaks=br )

## take other pi0, to match 


## let's try EM, with pi0 fixed at the value from fdrtool, but to estimate f_A on uniform
## scale and on bins determined in h

# initialize

fA <- h$density/length(h$density)  ### start at margin

notdone <- TRUE
imax <- 500
tol <- 1/10^4
nstep <- 1
nbin <- length(h$density) ## converting to masses
while(notdone)
 {
  fmarg <- (1-pi0)*fA + pi0*(1/nbin)
  p.tmp <- fA*(1-pi0)/fmarg
  p.update <- h$counts*p.tmp/sum( h$counts*p.tmp )
  notdone <- ! ( sum( abs( p.update-fA ) ) < tol | nstep > imax  )
  nstep <- nstep+1
  fA <- p.update
  #print(nstep)
 }

# that get's a bona fide estimate of distn given alternative, fixing the mixing rateo

## now map back to original stats

tg <- qt( h$breaks[2:(nbin+1)], df=22 )
pg <- fA*nbin  ## convert to density scale
tmp.fun <- approxfun(tg, pg )

### now I'm missing the Jacobian to get us back to density on the t-statistic scale
## simple trick
iii <- integrate(tmp.fun, min(tg), max(tg) )
FA.fun <- approxfun( tg, pg/iii$value )
 

## that seems to work [check plot( tt, FA.fun(tt), pch="." ), and the choice of bins leaves no NA's

lfdr <-  pi0 * f0.val/( pi0*f0.val + (1-pi0)*FA.fun(tt) )


## 
plot( lfdr, g$lfdr, pch="." )
abline(0,1, col="red" )
## there's a ton of information in the sidedness in this case
## lots of peptides with lfdr <= .05 by the t-stat method (constraining the pi0 and fA)
##  || sum(lfdr < .05) is 237; check > plot(tt, lfdr, pch="." )
##  || 28 less than .01
## compared to the method which used only p-values
## recall, Efron's locfdr would not produce any estimates using a fixed pi0 and 
## but it seems there is just a ton of increased (but small effect) antibodies, which makes
## sense owing to the selected proteins

fmarg <- pi0*dt( tg, df=22 ) + (1-pi0)*FA.fun(tg)
plot( tg, fmarg, lwd=2, xlab="t stat", ylab="constrained mixture", type="l" )
lines( tg, pi0*dt(tg, df=22 ) , col="green", lwd=2 )
## this shows that FA is a bad estimate, but mostly because pi0 is so badly estimated
## that no mixture will fit the margin well at that high pi0
## todo: rerun with a smaller pi0 just to show an improved marginal fit

## March 28
## now, it's not so clear why using two dimensions will get us a lot more than this 
## hand-crafted t-statistic approach.   Evidently we get stuff using mixtwice and ASH
## it would be constructive to compare that yeild with the yeild of this custom t approach
## using the same pi0 value


