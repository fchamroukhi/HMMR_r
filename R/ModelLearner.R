#' @export
emHMMR <- function(X, Y, K, p, variance_type = c("heteroskedastic", "homoskedastic"), n_tries = 1, max_iter = 1500, threshold = 1e-6, verbose = FALSE) {

  fData <- FData$new(X, Y)

  nb_good_try <- 0
  total_nb_try <- 0
  best_loglik <- -Inf
  cputime_total <- c()

  while (nb_good_try < n_tries) {
    start_time <- Sys.time()

    if (n_tries > 1) {
      print(paste("EM try n?", (nb_good_try + 1)))
    }
    total_nb_try <- total_nb_try + 1

    ## EM Initializaiton step
    ## Initialization of the Markov chain params, the regression coeffs, and the variance(s)
    variance_type <- match.arg(variance_type)
    param <- ParamHMMR$new(fData = fData, K = K, p = p, variance_type = variance_type)
    param$initHmmr(nb_good_try + 1)

    iter <- 0
    prev_loglik <- -Inf
    converged <- FALSE
    top <- 0

    stat <- StatHMMR(paramHMMR = param)

    while ((iter <= max_iter) && !converged) {

      ## E step : calculate tge tau_tk (p(Zt=k|y1...ym;theta)) and xi t_kl (and the log-likelihood) by
      #  forwards backwards (computes the alpha_tk et beta_tk)
      stat$EStep(param)

      ## M step
      param$MStep(stat)

      ## End of an EM iteration

      iter <-  iter + 1

      # test of convergence
      lambda <- 1e-9 # if a bayesian prior on the beta's
      stat$loglik <- stat$loglik + log(lambda)

      if (verbose) {
        print(paste('HMM_regression | EM   : Iteration :', iter, ' Log-likelihood : ', stat$loglik))
      }

      if ((prev_loglik - stat$loglik) > 1e-4) {
        top <- top + 1
        if (top == 10) {
          stop(print(paste('!!!!! The loglikelihood is decreasing from', prev_loglik, ' to ', stat$loglik)))
        }
      }

      converged <- (abs(stat$loglik - prev_loglik) / abs(prev_loglik) < threshold)
      if (is.na(converged)) {
        converged <- FALSE
      } # Basically for the first iteration when prev_loglik is Inf

      prev_loglik <- stat$loglik
      stat$stored_loglik[iter] <- stat$loglik

    } # END EM LOOP

    cputime_total[nb_good_try + 1] <- Sys.time() - start_time

    # at this point we have computed param and stat that contains all the information

    if (n_tries > 1) {
      print(paste('loglik_max = ', stat$loglik))
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
    print(paste('best_loglik:  ', statSolution$loglik))
  }

  # Smoothing state sequences : argmax(smoothing probs), and corresponding binary allocations partition
  statSolution$MAP()

  # FINISH computation of statSolution
  statSolution$computeStats(paramSolution, cputime_total)

  return(ModelHMMR(paramHMMR = paramSolution, statHMMR = statSolution))
}
