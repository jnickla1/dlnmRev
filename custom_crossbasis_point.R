custom_crossbasis_point <- function (varpoints,lagpoints,daypoints, argvar = list(), arglag = list()) 
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

  basislag2 <- do.call("onebasis", modifyList(arglag, list(x = lagpoints)))
  basisday2 <- do.call("onebasis", modifyList(arglag, list(x = daypoints)))
  basisvar2 <- do.call("onebasis", modifyList(argvar, list(x = varpoints)))
  npoints = length(varpoints)
  nbasis =  dim(basisvar2)[2]*dim(basislag2)[2]*dim(basislag2)[2]
  crossbasisreturnM <- array(rep(0, nbasis*npoints), c(npoints,nbasis))

  for ( point in seq(length = npoints)) {
    
    crossbasisreturnM[point,] =  basisvar2[point,] %x% basislag2[point,] %x% basisday2[point,] 
  }
  
  return(crossbasisreturnM)
}
