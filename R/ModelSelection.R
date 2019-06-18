#' @export
selectHMMR <- function(X, Y, Kmin = 1, Kmax = 10, pmin = 0, pmax = 4, criterion = c("BIC", "AIC")) {

  criterion <- match.arg(criterion)

  vhmmr <- Vectorize(function(K, p, X1 = X, Y1 = Y) emHMMR(X = X1, Y = Y1, K, p),
                     vectorize.args = c("K", "p"))

  hmmr <- outer(Kmin:Kmax, pmin:pmax, vhmmr)

  if (criterion == "BIC") {
    results <- apply(hmmr, 1:2, function(x) x[[1]]$statHMMR$BIC)
  } else {
    results <- apply(hmmr, 1:2, function(x) x[[1]]$statHMMR$AIC)
  }
  rownames(results) <- sapply(Kmin:Kmax, function(x) paste0("(K = ", x, ")"))
  colnames(results) <- sapply(pmin:pmax, function(x) paste0("(p = ", x, ")"))


  selected <- hmmr[which(results == max(results), arr.ind = T)][[1]]

  cat(paste0("The HMMR model selected via the \"", criterion, "\" has K = ",
             selected$paramHMMR$K, " regimes \n and the order of the ",
             "polynomial regression is p = ", selected$paramHMMR$p, "."))
  cat("\n")
  cat(paste0("BIC = ", selected$statHMMR$BIC, "\n"))
  cat(paste0("AIC = ", selected$statHMMR$AIC, "\n"))

  return(selected)

}
