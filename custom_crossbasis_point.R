custom_crossbasis_point <- function (varpoints,lagpoints,daypoints,xminmax, lag, argvar = list(), arglag = list()) 
{ #generate a 1d basis array, where the first dimension is the identity of each basis function
  #checkcrossbasis(argvar, arglag, list(...))
  
  #argvar0 <- list(fun="strata",breaks=varbreaks, intercept=TRUE) # the groups of counts in the individual x's along var
  
  errormes <- "arguments 'varpoints, lagpoints, daypoints' not consistent."
  if (is.vector(varpoints) && is.vector(lagpoints)&& is.vector(daypoints)) {
    if ( length(varpoints) != length(lagpoints)) 
      stop(errormes)
    else if ( length(varpoints) != length(daypoints)) 
      stop(errormes)}
  else
  stop(errormes)

  basislag2a <- do.call("onebasis", modifyList(arglag, list(x = c(0,lag,lagpoints))))
  basisday2a <- do.call("onebasis", modifyList(arglag, list(x = c(0,lag,daypoints))))
  basisvar2a <- do.call("onebasis", modifyList(argvar, list(x = c(xminmax[1:2],varpoints))))
  basislag2 <- basislag2a[c(-1,-2),]
  basisday2 <- basisday2a[c(-1,-2),]
  basisvar2 <- basisvar2a[c(-1,-2),]
  npoints = length(varpoints)
  nbasis =  dim(basisvar2)[2]*dim(basislag2)[2]*dim(basislag2)[2]
  crossbasisreturnM <- array(rep(0, nbasis*npoints), c(npoints,nbasis))

  for ( point in seq(length = npoints)) {
    
    crossbasisreturnM[point,] =  basisvar2[point,] %x% basislag2[point,] %x% basisday2[point,] 
  }
  
  return(crossbasisreturnM)
}
