lognorm_params <- function(xmean, xmin, xmax, percentile=0.05){
  logins = log(c(xmean, xmin, xmax))
  meanlog = mean(logins)
  widesdlog = (logins[3] - logins[2])/2
  sdlog= widesdlog/-qnorm(percentile/2,0,1)
  return(c(meanlog,sdlog))
}

lognorm_params_kelly <- function(xmean, x10, x90, percentile=0.1){
  logins = log(c(xmean, x10, x90))
  meanlog0 = logins[1]
  widesdlog = (logins[3] - logins[2])/2
  sdlog0= widesdlog/-qnorm(percentile,0,1)
  alpha = logins[3] + logins[2] - 2*logins[1]
  delta = alpha / sqrt(1+alpha^2)
  omega = sdlog0 / sqrt(1-2*delta^2/pi)
  eta = meanlog0 -omega * delta * sqrt(2/pi)
  return(c(eta,omega,alpha))
}

PmeasB_below_LNassumed <- function(exceeds,meanlg,sdlg, alph, beta,  P){
  f <- function(y) dbeta(y, alph, beta)*(1- plnorm(y+exceeds, meanlog = meanlg, sdlog = sdlg ))
  integ = integrate(f, 0,1-exceeds)
  return( integ$value - P)
}

PmeasB_times_LNassumed <- function(exceeds,meanlg,sdlg, alph, beta,  P){
  f <- function(y) dbeta(y, alph, beta)*(1- plnorm(y*exceeds, meanlog = meanlg, sdlog = sdlg ))
  integ = integrate(f, 0,1)
  return( integ$value - P)
}

#uniroot(PmeasB_below_LNassumed , c(-1,1), extendInt="downX")
#unr <- uniroot(PmeasB_below_LNassumed , c(-1,1), extendInt="downX",meanlg=-2,sdlg=.1, alph=2, beta=18,  P=0.5)
#log(unr$root)


library("sn")#skewnormal
plsnorm <- function( x, eta, omega, alpha){
  return(pnorm((log(x)-eta)/omega) - 2*T.Owen((log(x)-eta)/omega,alpha))
}

PmeasB_times_LNskew <- function(exceeds,meanlg,sdlg,skew, alph, beta,  P){
  f <- function(y) dbeta(y, alph, beta)*(1- plsnorm(x = y*exceeds, eta = meanlg, omega = sdlg ,alpha=skew))
  integ = integrate(f, 0,1,rel.tol = 1e-8)
  return( integ$value - P)
}

CmeasB_times_LNskewfit <- function(percdata, alph, beta,  P=0.5){
  parms <- lognorm_params_kelly(percdata[1],percdata[2],percdata[3])
  unr <- uniroot(PmeasB_times_LNskew, c(0.0001,500), extendInt="downX",meanlg=parms[1],sdlg=parms[2],skew=parms[3], alph, beta,P)
  return(-log(unr$root))
}

