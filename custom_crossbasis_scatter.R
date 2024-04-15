custom_crossbasis_scatter <- function (x, lag, varbreaks, group = NULL, pslags = NULL, dayweights = NULL,densReturn=FALSE) 
{
  #checkcrossbasis(argvar, arglag, list(...))
  #argvar <- list(fun="strata",breaks=varbreaks, intercept=TRUE)
  if(is.null(dayweights)){
    dayweights <-  matrix(1,length(x))}
    
  lag <- if (missing(lag)) c(0, NCOL(x) - 1) else mklag(lag)
  x <- as.matrix(x)
  dim <- dim(x)
  if (!dim[2] %in% c(1L, diff(lag) + 1L)) 
    stop("NCOL(x) must be equal to 1 (if x is a time series vector), ", 
         "otherwise to the lag period (for x as a matrix of lagged occurrences)")
  #basisvar <- do.call("onebasis", modifyList(argvar, list(x = as.numeric(x)))) #hate this weird implementation
  basisvar <- matrix(0,length(x),length(varbreaks)+1) #reduce the dimension of the variable x to argvar dimension
  onelocs <- findInterval(x, varbreaks)
  for (i in seq(length(x))){ basisvar[i,onelocs[i]+1] <- 1}
  
  #basislag <- do.call("onebasis", modifyList(arglag, list(x = seqlag(lag))))
  
  #if (!is.null(group))
  #  checkgroup(group, x, basisvar, lag)
  nlag_days = lag[2]+1
  crossbasisscatter <- array(rep(0, dim(basisvar)[2]*nlag_days*nlag_days), c(dim(basisvar)[2],nlag_days,nlag_days))
  
  fines = seq(min(x),max(x),0.1)
  #fines_b = do.call("onebasis", modifyList(argvar, list(x = fines)))
  fines_b <- matrix(0,length(fines),length(varbreaks)+1) #reduce the dimension of the variable x to argvar dimension
  onelocsf <- findInterval(fines, varbreaks)
  for (i in seq(length(fines))){ fines_b[i,onelocsf[i]+1] <- 1}
  
  if (densReturn){
    norms_density = colSums(fines_b)*0.1 #density per 0.1 degrees
  }else{
    norms_density = matrix(1,length(varbreaks)+1) #count number of observations in each bin
  }
  
  mat <- matrix(data=NA,nrow=dim[1],ncol=lag[2]+1)
  colnames(mat) <- seq(0,lag[2])
  gs <- table(group)
  for (v in seq(length = ncol(basisvar))) {
    if (dim[2] == 1L) { #data is just a list of values with group identifying which belongs to each patient
      st<-1L 
      for (gi in gs){
      en<- st+ gi -1 #last entry for each patient, st is the first
      mat[st:en, 1:gi ] <- as.matrix(Lag(basisvar[st:en, v] , seqlag(c(0, gi - 1)))) *dayweights[st:en]
      #store each lag day in the past for each patient, with columns being these lag days.
      st<- en + 1L
      }
    }
    else mat <- matrix(basisvar[, v], ncol = diff(lag) + 1)
    mat[is.na(mat)] <- 0 #zeros where nothing was reached yet
    #for (d in seq(length = ncol(basislag))) {
    #  dbasislag <- basislag[as.integer(pslags), d] #compute lags
   # for (l in seq(length = ncol(basislag))) {
      #compute the influence of each crossbasis function on the computed history*lag triangle for each patient
   #   crossbasis[, ncol(basislag)*(ncol(basislag) * (v - 1) + l-1) + d] <- dbasislag * (mat %*% (basislag[, l]))
   # }}
    
    for (r in seq_len(dim[1])){
      rlag=pslags[r]
      crossbasisscatter[v,rlag,]  = crossbasisscatter[v,rlag,] + mat[r,] / max(norms_density[v],1)
    }

  }

 # ind <- match(names(formals(attributes(basisvar)$fun)), names(attributes(basisvar)), nomatch = 0)
 # argvar <- c(attributes(basisvar)["fun"], attributes(basisvar)[ind])

 # argvar$cen <- attributes(basisvar)$cen
 # attributes(crossbasisscatter) <- c(attributes(crossbasis), list(df = c(ncol(basisvar)), 
 #                               range = range(x, na.rm = T), lag = lag, argvar = argvar, arglag = arglag))
  #if (!is.null(group)) 
 #   attributes(crossbasisscatter)$group <- length(unique(group))
  #class(crossbasis) <- c("crossbasis", "matrix")

  
  return(crossbasisscatter)
}
#library('plot.matrix')
#plot(totals[2,,],breaks=seq(700,800,20))