#' emHMMR is used to fit a HMMR model.
#'
#' emHMMR is used to fit a HMMR model. The estimation method is performed by
#' the Expectation-Maximization algorithm.
#'
#' @details emHMMR function is based on the EM algorithm. This function starts
#' with an initialization of the parameters done by the method `initParam` of
#' the class [ParamHMMR][ParamHMMR], then it alternates between a E-Step
#' (method of the class [StatHMMR][StatHMMR]) and a M-Step (method of the class
#' [ParamHMMR][ParamHMMR]) until convergence (until the absolute difference of
#' log-likelihood between two steps of the EM algorithm is less than the
#' `threshold` parameter).
#'
#' @param X Numeric vector of length \emph{m} representing the covariates/inputs
#' \eqn{x_{1},\dots,x_{m}}.
#' @param Y Numeric vector of length \emph{m} representing the observed
#' response/output \eqn{y_{1},\dots,y_{m}}.
#' @param K The number of regimes (mixture components).
#' @param p The order of the polynomial regression.
#' @param variance_type Optional character indicating if the model is
#' "homoskedastic" or "heteroskedastic". By default the model is
#' "heteroskedastic".
#' @param n_tries Optional. Number of times EM algorithm will be launched.
#' The solution providing the highest log-likelihood will be returned.
#'
#' If `n_tries` > 1, then for the first pass, parameters are initialized
#' by uniformly segmenting the data into K segments, and for the next passes,
#' parameters are initialized by randomly segmenting the data into K contiguous
#'  segments.
#' @param max_iter Optional. The maximum number of iterations for the EM algorithm.
#' @param threshold Optional. A numeric value specifying the threshold for the relative
#'  difference of log-likelihood between two steps  of the EM as stopping
#'  criteria.
#' @param verbose Optional. A logical value indicating whether values of the
#' log-likelihood should be printed during EM iterations.
#' @return EM returns an object of class [ModelHMMR][ModelHMMR].
#' @seealso [ModelHMMR], [ParamHMMR], [StatHMMR]
#' @export
emHMMR <- function(X, Y, K, p = 3, variance_type = c("heteroskedastic", "homoskedastic"), n_tries = 1, max_iter = 1500, threshold = 1e-6, verbose = FALSE) {

  nb_good_try <- 0
  total_nb_try <- 0
  best_loglik <- -Inf
  cputime_total <- c()

  while (nb_good_try < n_tries) {
    start_time <- Sys.time()

    if (n_tries > 1) {
      if (verbose) {
        cat(paste0("EM try number: ", nb_good_try + 1, "\n\n"))
      }
    }
    total_nb_try <- total_nb_try + 1

    ## EM Initializaiton step
    ## Initialization of the Markov chain params, the regression coeffs, and the variance(s)
    variance_type <- match.arg(variance_type)
    param <- ParamHMMR$new(X = X, Y = Y, K = K, p = p, variance_type = variance_type)
    param$initParam(nb_good_try + 1)

    iter <- 0
    prev_loglik <- -Inf
    converged <- FALSE
    top <- 0

    stat <- StatHMMR(paramHMMR = param)

    while ((iter <= max_iter) && !converged) {

      # E step : calculate tge tau_tk (p(Zt=k|y1...ym;theta)) and xi t_kl (and the log-likelihood) by
      #  forwards backwards (computes the alpha_tk et beta_tk)
      stat$EStep(param)

      # M step
      param$MStep(stat)

      # End of an EM iteration

      iter <-  iter + 1

      # Test of convergence
      lambda <- 1e-9 # If a bayesian prior on the beta's
      stat$loglik <- stat$loglik + log(lambda)

      if (verbose) {
        cat(paste0("EM: Iteration : ", iter, " || log-likelihood : "  , stat$loglik, "\n"))
      }

      if ((prev_loglik - stat$loglik) > 1e-4) {
        top <- top + 1
        if (top == 10) {
          warning(paste0("EM log-likelihood is decreasing from ", prev_loglik, "to ", stat$loglik, " !"))
        }
      }

      converged <- (abs(stat$loglik - prev_loglik) / abs(prev_loglik) < threshold)
      if (is.na(converged)) {
        converged <- FALSE
      } # Basically for the first iteration when prev_loglik is Inf

      prev_loglik <- stat$loglik
      stat$stored_loglik[iter] <- stat$loglik

    } # End of the EM loop

    cputime_total[nb_good_try + 1] <- Sys.time() - start_time

    if (n_tries > 1) {
      cat(paste0("Max value of the log-likelihood: ", stat$loglik, "\n"))
    }

    if (length(param$beta) != 0) {
      nb_good_try <- nb_good_try + 1
      total_nb_try <- 0

      if (stat$loglik > best_loglik) {

        statSolution <- stat$copy()
        paramSolution <- param$copy()

        best_loglik <- stat$loglik
      }
    }

    if (total_nb_try > 500) {
      stop(paste("can't obtain the requested number of classes"))
    }

  }

  if (n_tries > 1) {
    cat(paste0("Best value of the log-likelihood: ", statSolution$loglik, "\n"))
  }

  # Smoothing state sequences : argmax(smoothing probs), and corresponding binary allocations partition
  statSolution$MAP()

  # Finish the computation of statistics
  statSolution$computeStats(paramSolution, cputime_total)

  return(ModelHMMR(paramHMMR = paramSolution, statHMMR = statSolution))
}
