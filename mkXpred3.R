mkXpred3 <- function (type, basis, at, predvar, predlag, cen) 
{
  varvec <- rep(at, length(predlag)*length(predlag)) #note this grouping is flipped to allow matrix creation 
  lagvec <- rep(predlag, each = length(predvar), times = length(predlag))
  dagvec <- rep(predlag, each = length(predvar)*length(predlag))
  if (type == "cb") {
    basisvar <- do.call("onebasis", c(list(x = varvec), 
                                      attr(basis, "argvar")))
    basislag <- do.call("onebasis", c(list(x = lagvec), 
                                      attr(basis, "arglag")))
    basislagd <- do.call("onebasis", c(list(x = dagvec), 
                                      attr(basis, "arglag")))
    if (!is.null(cen)) {
      basiscen <- do.call("onebasis", c(list(x = cen), 
                                        attr(basis, "argvar")))
      basisvar <- scale(basisvar, center = basiscen, scale = FALSE)
    }
    Xpred <- tensor.prod.model.matrix(list(basisvar, basislag, basislagd))
  }
  return(Xpred)
}