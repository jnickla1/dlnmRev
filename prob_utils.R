lognorm_params <- function(xmean, xmin, xmax, percentile=0.05){
  logins = log(c(xmean, xmin, xmax))
  meanlog = mean(logins)
  widesdlog = (logins[3] - logins[2])/2
  sdlog= widesdlog/-qnorm(percentile/2,0,1)
  return(c(meanlog,sdlog))
}

PmeasB_below_LNassumed <- function(exceeds,meanlg,sdlg, alph, beta,  P){
  f <- function(y) dbeta(y, alph, beta)*(1- plnorm(y+exceeds, meanlog = meanlg, sdlog = sdlg ))
  integ = integrate(f, 0,1-exceeds)
  return( integ$value - P)
}

#uniroot(PmeasLN_exceeds_bassumed, c(-1,1), extendInt="upX")