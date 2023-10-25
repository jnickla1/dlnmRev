custom_crosspred <- function (basis, model = NULL, coef = NULL, vcov = NULL, model.link = NULL, 
    at = NULL, from = NULL, to = NULL, by = NULL, lag, bylag = 1, 
    cen = NULL, ci.level = 0.95, cumul = FALSE) 
{   source("~/Documents/dlnmRev/mkXpred3.R")
    environment(mkXpred3) <- asNamespace('dlnm')
    type <- if (any(class(basis) %in% "crossbasis")) 
        "cb"
    else if (any(class(basis) %in% "onebasis")) 
        "one"
    else "gam"
    errormes <- "arguments 'basis' and 'model' not consistent. See help(crosspred)"
    if (type == "gam") {
        if (!is.character(basis) || length(basis) > 1L) 
            stop(errormes)
        if (is.null(model) || !any(class(model) %in% "gam")) 
            stop(errormes)
        name <- basis
        sterms <- sapply(model$smooth, function(x) x$term[1])
        if (name %in% sterms) 
            basis <- model$smooth[[which(sterms == name)[1]]]
        else stop(errormes)
        if (length(which(sterms == name)) > 1) 
            warning(paste(name, "included in multiple smoothers, only the first one taken"))
        if (!"cb.smooth" %in% class(basis) && basis$dim > 1L) 
            stop("predictions not provided for multi-dimensional smoothers other than 'cb'")
    }
    else name <- deparse(substitute(basis))
    origlag <- switch(type, cb = attr(basis, "lag"), one = c(0, 
        0), gam = if (is.null(basis$lag)) c(0, 0) else basis$lag)
    lag <- if (missing(lag)) 
        origlag
    else mklag(lag)
    if (!all(lag == origlag) && cumul) 
        stop("cumulative prediction not allowed for lag sub-period")
    lagfun <- switch(type, cb = attr(basis, "arglag")$fun, one = NULL, 
        gam = if (basis$dim == 1L) NULL else basis$margin[[2]]$fun)
    if (bylag != 1L && !is.null(lagfun) && lagfun == "integer") 
        stop("prediction for non-integer lags not allowed for type 'integer'")
    if (is.null(model) && (is.null(coef) || is.null(vcov))) 
        stop("At least 'model' or 'coef'-'vcov' must be provided")
    if (!is.numeric(ci.level) || ci.level >= 1 || ci.level <= 
        0) 
        stop("'ci.level' must be numeric and between 0 and 1")
    cond <- if (type == "gam") 
        with(basis, first.para:last.para)
    else if (ncol(basis) == 1L) 
        name
    else if (type == "one") 
        paste(name, "[[:print:]]*b[0-9]{1,2}", sep = "")
    else paste(name, "[[:print:]]*v[0-9]{1,2}\\.l[0-9]{1,2}", 
        sep = "")
    if (!is.null(model)) {
        model.class <- class(model)
        coef <- getcoef(model, model.class)
        vcov <- getvcov(model, model.class)
        indcoef <- if (type == "gam") 
            cond
        else grep(cond, names(coef))
        indvcov <- if (type == "gam") 
            cond
        else grep(cond, rownames(vcov))
        coef <- coef[indcoef]
        vcov <- vcov[indvcov, indvcov, drop = FALSE]
        model.link <- getlink(model, model.class, model.link)
    }
    else model.class <- NA
    npar <- if (type == "gam") 
        length(indcoef)
    else ncol(basis)
    coef[is.na(coef)] <- 0
    vcov[is.na(vcov)] <- 0
    if (length(coef) != npar || length(coef) != dim(vcov)[1] )# || need to deal with the nan coeff values at some point
       # any(is.na(coef)) || any(is.na(vcov))) 
        stop("coef/vcov not consistent with basis matrix. See help(crosspred)")
    range <- if (type == "gam") 
        range(model$model[[basis$term[1]]])
    else attr(basis, "range")
    at <- mkat(at, from, to, by, range, lag, bylag)
    predvar <- if (is.matrix(at)) 
        rownames(at)
    else at
    predlag <- seqlag(lag, bylag)
    cen <- mkcen(cen, type, basis, range) #something funny happening here, not calculating correctly
    if (type == "one") 
        attributes(basis)$cen <- NULL
    if (type == "cb") 
        attributes(basis)$argvar$cen <- NULL
    Xpred <- mkXpred3(type, basis, at, predvar, predlag, cen)
    matfit <- matrix(Xpred %*% coef, length(predvar), length(predlag)*length(predlag))
    matse <- matrix(sqrt(pmax(0, rowSums((Xpred %*% vcov) * Xpred))), 
        length(predvar), length(predlag)*length(predlag))
    rownames(matfit) <- rownames(matse) <- predvar
    colnames(matfit) <- colnames(matse) <- outer("lag", predlag, 
        paste, sep = "")
    predlag <- seqlag(lag)
    Xpred <- mkXpred3(type, basis, at, predvar, predlag, cen)
    Xpredall <- 0
    if (cumul) {
        cumfit <- cumse <- matrix(0, length(predvar), length(predlag))
    }
    for (i in seq(length(predlag))) {
        ind <- seq(length(predvar)) + length(predvar) * (i - 
            1)
        Xpredall <- Xpredall + Xpred[ind, , drop = FALSE]
        if (cumul) {
            cumfit[, i] <- Xpredall %*% coef
            cumse[, i] <- sqrt(pmax(0, rowSums((Xpredall %*% 
                vcov) * Xpredall)))
        }
    }
    allfit <- as.vector(Xpredall %*% coef)
    allse <- sqrt(pmax(0, rowSums((Xpredall %*% vcov) * Xpredall)))
    names(allfit) <- names(allse) <- predvar
    if (cumul) {
        rownames(cumfit) <- rownames(cumse) <- predvar
        colnames(cumfit) <- colnames(cumse) <- outer("lag", seqlag(lag), 
            paste, sep = "")
    }
    list <- list(predvar = predvar)
    if (!is.null(cen)) 
        list$cen <- cen
    list <- c(list, list(lag = lag, bylag = bylag, coefficients = coef, 
        vcov = vcov, matfit = matfit, matse = matse, allfit = allfit, 
        allse = allse))
    if (cumul) 
        list <- c(list, list(cumfit = cumfit, cumse = cumse))
    z <- qnorm(1 - (1 - ci.level)/2)
    if (!is.null(model.link) && model.link %in% c("log", "logit")) {
        list$matRRfit <- exp(matfit)
        list$matRRlow <- exp(matfit - z * matse)
        list$matRRhigh <- exp(matfit + z * matse)
        list$allRRfit <- exp(allfit)
        list$allRRlow <- exp(allfit - z * allse)
        names(list$allRRlow) <- names(allfit)
        list$allRRhigh <- exp(allfit + z * allse)
        names(list$allRRhigh) <- names(allfit)
        if (cumul) {
            list$cumRRfit <- exp(cumfit)
            list$cumRRlow <- exp(cumfit - z * cumse)
            list$cumRRhigh <- exp(cumfit + z * cumse)
        }
    }
    else {
        list$matlow <- matfit - z * matse
        list$mathigh <- matfit + z * matse
        list$alllow <- allfit - z * allse
        names(list$alllow) <- names(allfit)
        list$allhigh <- allfit + z * allse
        names(list$allhigh) <- names(allfit)
        if (cumul) {
            list$cumlow <- cumfit - z * cumse
            list$cumhigh <- cumfit + z * cumse
        }
    }
    list$ci.level <- ci.level
    list$model.class <- model.class
    list$model.link <- model.link
    class(list) <- "crosspred"
    return(list)
}