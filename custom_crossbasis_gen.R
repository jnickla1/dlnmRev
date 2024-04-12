custom_crossbasis_gen <- function (xminmax, lag,varbreaks, argvar = list(), arglag = list()) 
{ #generate a 4d basis array, where the first dimension is the identity of each basis function
  #checkcrossbasis(argvar, arglag, list(...))
  
  #argvar0 <- list(fun="strata",breaks=varbreaks, intercept=TRUE) # the groups of counts in the individual x's along var
  x=seq(xminmax[1],xminmax[2],xminmax[3])
  
  lag <- if (missing(lag)) 
    c(0, NCOL(x) - 1)
  else mklag(lag)
  
  #basisvar <- do.call("onebasis", modifyList(argvar, list(x = as.numeric(x)))) #hate this weird implementation
  basisvar <- matrix(0,length(x),length(varbreaks)+1) #reduce the dimension of the variable x to argvar dimension
  onelocs <- findInterval(x, varbreaks)
  for (i in seq(length(x))){ basisvar[i,onelocs[i]+1] <- 1}
  
  if (length(arglag) == 0L || diff(lag) == 0L) 
    arglag <- list(fun = "strata", df = 1, intercept = TRUE)
  if ((is.null(arglag$fun) || "intercept" %in% names(formals(arglag$fun))) && 
      sum(pmatch(names(arglag), "intercept", nomatch = 0)) == 0) 
    arglag$intercept <- TRUE
  arglag$cen <- NULL
  basislag2 <- do.call("onebasis", modifyList(arglag, list(x = seqlag(lag)))) #new possibly spline basis from all-integer lag sequence
  basisvar2 <- do.call("onebasis", modifyList(argvar, list(x = as.numeric(x))))
  
  varbasistransfer <- matrix(data=NA,nrow=dim(basisvar)[2],ncol=dim(basisvar2)[2])
  
  #if (!is.null(group))
  #  checkgroup(group, x, basisvar, lag)
  nlag_days = lag[2]+1
  nbasis =  dim(basisvar2)[2]*dim(basislag2)[2]*dim(basislag2)[2]
  crossbasisdesignM <- array(rep(0, nbasis*dim(basisvar)[2]*nlag_days*nlag_days), c(nbasis,dim(basisvar)[2],nlag_days,nlag_days))
  #has later dimensions: origbasisvar, lag, day

  for ( bvar2 in seq(length = dim(basisvar2)[2])) {
    origbasis_weight = t(basisvar2[,bvar2]) %*% basisvar / t(rep(1, dim(basisvar)[1])) %*% basisvar 
    #average impact from that col of basisvar2 on each basis of basisvar
    for ( blag2 in seq(length = dim(basislag2)[2])) {
      for ( bday2 in seq(length = dim(basislag2)[2])) {
        v  =  (bvar2-1) * (dim(basislag2)[2]*dim(basislag2)[2]) +  (blag2 -1) * dim(basislag2)[2] + bday2
        for (c in seq(dim(origbasis_weight)[2])){
        crossbasisdesignM[v,c,,] = origbasis_weight[1,c]* array(kronecker(basislag2[,blag2],basislag2[,bday2]),c(nlag_days,nlag_days))
        }
          #kronecker(array(kronecker(origbasis_weight,basislag2[,blag2]),dim=c(dim(basisvar)[2],nlag_days,1)),basislag2[,bday2])
        
      }}}
  
  for ( d in seq(2,nlag_days)){
    crossbasisdesignM[,,0:(d-1),d] = 0 
  }
    #for (d in seq(length = ncol(basislag))) {
    #  dbasislag <- basislag[as.integer(pslags), d] #compute lags
   # for (l in seq(length = ncol(basislag))) {
      #compute the influence of each crossbasis function on the computed history*lag triangle for each patient
   #   crossbasis[, ncol(basislag)*(ncol(basislag) * (v - 1) + l-1) + d] <- dbasislag * (mat %*% (basislag[, l]))
   # }}

 # ind <- match(names(formals(attributes(basisvar)$fun)), names(attributes(basisvar)), nomatch = 0)
 # argvar <- c(attributes(basisvar)["fun"], attributes(basisvar)[ind])

 # argvar$cen <- attributes(basisvar)$cen
 # attributes(crossbasisscatter) <- c(attributes(crossbasis), list(df = c(ncol(basisvar)), 
 #                               range = range(x, na.rm = T), lag = lag, argvar = argvar, arglag = arglag))
  #if (!is.null(group)) 
 #   attributes(crossbasisscatter)$group <- length(unique(group))
  #class(crossbasis) <- c("crossbasis", "matrix")
  return(crossbasisdesignM)
}
#library('plot.matrix')
#plot(totals[2,,],breaks=seq(700,800,20))