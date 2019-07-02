#' selectHMMR implements the model selection procedure.
#'
#' @details selectHMMR implements the model selection. This function runs every
#'   HMMR model by varying the number of regimes `K` from `Kmin` to `Kmax` and
#'   the order of the polynomial regression `p` from `pmin` to `pmax`. The model
#'   having the highest value of the chosen selection criterion is then
#'   selected.
#'
#' @param X Numeric vector of length \emph{m} representing the covariates/inputs
#'   \eqn{x_{1},\dots,x_{m}}.
#' @param Y Numeric vector of length \emph{m} representing the observed
#'   response/output \eqn{y_{1},\dots,y_{m}}.
#' @param Kmin The minimum number of regimes (c components).
#' @param Kmax The maximum number of regimes (HMMR components).
#' @param pmin The minimum order of the polynomial regression.
#' @param pmax The maximum order of the polynomial regression.
#' @param criterion The criterion used to select the HMMR model ("BIC", "AIC").
#' @return selectHMMR returns an object of class [ModelHMMR][ModelHMMR]
#'   representing the selected HMMR model according to the chosen `criterion`.
#' @seealso [ModelHMMR]
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
