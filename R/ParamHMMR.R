#' @export
ParamHMMR <- setRefClass(
  "ParamHMMR",
  fields = list(
    fData = "FData",
    phi = "matrix",

    K = "numeric", # Number of regimes
    p = "numeric", # Dimension of beta (order of polynomial regression)
    variance_type = "numeric",
    nu = "numeric", # Degree of freedom

    prior = "matrix",
    trans_mat = "matrix",
    beta = "matrix",
    sigma = "matrix",
    mask = "matrix"
  ),
  methods = list(

    initialize = function(fData = FData(numeric(1), matrix(1)), K = 2, p = 2, variance_type = 1) {

      fData <<- fData

      phi <<- designmatrix(x = fData$X, p = p)$XBeta

      K <<- K
      p <<- p
      variance_type <<- variance_type

      if (variance_type == variance_types$homoskedastic) {
        nu <<- K - 1 + K * (K - 1) + K * (p + 1) + 1
      }
      else{
        nu <<- K - 1 + K * (K - 1) + K * (p + 1) + K
      }

      prior <<- matrix(NA, ncol = K - 1)
      trans_mat <<- matrix(NA, K, K)
      beta <<- matrix(NA, p + 1, K)
      if (variance_type == variance_types$homoskedastic) {
        sigma <<- matrix(NA)
      }
      else{
        sigma <<- matrix(NA, K)
      }
      mask <<- matrix(NA, K, K)

    },

    initHmmr = function(try_algo = 1) {
      # function hmmr =  initHmmr(X, y, K, type_variance, EM_try)
      # initHmmr initialize the parameters of a Gaussian Hidden Markov Model
      # Regression (HMMR) model
      #
      # Inputs :
      #
      #           X: [nx(p+1)] regression desing matrix
      #           K : Number of polynomial regression components (regimes)
      #          	type_variance: hoskedastoc or heteroskedastic
      #           EM_try: number of the current EM run
      #
      # Outputs :
      #
      #         hmmr: the initial HMMR model. a structure composed of:
      #
      #         prior: [Kx1]: prior(k) = Pr(z_1=k), k=1...K
      #         trans_mat: [KxK], trans_mat(\ell,k) = Pr(z_t = k|z_{t-1}=\ell)
      #         reg_param: the paramters of the regressors:
      #                 betak: regression coefficients
      #                 sigma2k (or sigma2) : the variance(s). sigma2k(k) = variance of y(t) given z(t)=k; sigma2k(k) =
      #         sigma^2_k.
      #         and some stats: like the mask for a segmental model
      #
      #
      ##################################################################################

      # Initialization taking into account the constraint:

      # Initialization of the transition matrix
      maskM <- 0.5 * diag(K) # mask of order 1

      for (k in 1:(K - 1)) {
        ind <- which(maskM[k, ] != 0)
        maskM[k, ind + 1] <- 0.5
      }
      trans_mat <<- maskM
      mask <<- maskM

      # Initialization of the initial distribution
      prior <<- matrix(c(1, rep(0, K - 1)))

      # Initialization of regression coefficients and variances
      initHmmrRegressors(try_algo)

    },

    initHmmrRegressors = function(try_algo = 1) {

      if (try_algo == 1) { # Uniform segmentation into K contiguous segments, and then a regression

        zi <- round(fData$m / K) - 1

        s <- 0 # If homoskedastic
        for (k in 1:K) {
          yk <- fData$Y[((k - 1) * zi + 1):(k * zi)]
          Xk <- as.matrix(phi[((k - 1) * zi + 1):(k * zi), ])

          beta[, k] <<- solve(t(Xk) %*% Xk + (10 ^ -4) * diag(p + 1)) %*% t(Xk) %*% yk

          muk <- Xk %*% beta[, k]
          sk <- t(yk - muk) %*% (yk - muk)
          if (variance_type == variance_types$homoskedastic) {
            s <- (s + sk)
            sigma <<- s / fData$m
          }
          else {
            sigma[k] <<- sk / length(yk)
          }
        }
      }
      else{# Random segmentation into contiguous segments, and then a regression

        Lmin <- p + 1 + 1 # Minimum length of a segment
        tk_init <- rep(0, K)
        tk_init <- t(tk_init)
        tk_init[1] <- 0
        K_1 <- K
        for (k in 2:K) {
          K_1 <- K_1 - 1
          temp <- seq(tk_init[k - 1] + Lmin, fData$m - K_1 * Lmin)
          ind <- sample(1:length(temp), length(temp))
          tk_init[k] <- temp[ind[1]]
        }
        tk_init[K + 1] <- fData$m

        s <- 0
        for (k in 1:K) {
          i <- tk_init[k] + 1
          j <- tk_init[k + 1]
          yk <- fData$Y[i:j]
          Xk <- phi[i:j, ]
          beta[, k] <<- solve(t(Xk) %*% Xk + 1e-4 * diag(p + 1)) %*% t(Xk) %*% yk
          muk <- Xk %*% beta[, k]
          sk <- t(yk - muk) %*% (yk - muk)

          if (variance_type == variance_types$homoskedastic) {
            s <- s + sk
            sigma <<- s / fData$m

          }
          else{
            sigma[k] <<- sk / length(yk)
          }
        }
      }
    },

    MStep = function(statHMMR) {
      # Updates of the Markov chain parameters
      # Initial states prob: P(Z_1 = k)
      prior <<- matrix(normalize(statHMMR$tau_tk[1, ])$M)

      # Transition matrix: P(Zt=i|Zt-1=j) (A_{k\ell})
      trans_mat <<- mkStochastic(apply(statHMMR$xi_tkl, c(1, 2), sum))

      # For segmental HMMR: p(z_t = k| z_{t-1} = \ell) = zero if k<\ell (no back) of if k >= \ell+2 (no jumps)
      trans_mat <<- mkStochastic(mask * trans_mat)
      # Update of the regressors (reg coefficients betak and the variance(s) sigma2k)

      s <- 0 # If homoskedastic
      for (k in 1:K) {
        weights <- statHMMR$tau_tk[, k]

        nk <- sum(weights) # Expected cardinal number of state k
        Xk <- phi * (sqrt(weights) %*% matrix(1, 1, p + 1)) # [n*(p+1)]
        yk <- fData$Y * sqrt(weights) # dimension: [(nx1).*(nx1)] = [nx1]

        # Regression coefficients
        lambda <- 1e-9 # If a bayesian prior on the beta's
        bk <- (solve(t(Xk) %*% Xk + lambda * diag(p + 1)) %*% t(Xk)) %*% yk
        beta[, k] <<- bk

        # Variance(s)
        z <- sqrt(weights) * (fData$Y - phi %*% bk)
        sk <- t(z) %*% z
        if (variance_type == variance_types$homoskedastic) {
          s <- (s + sk)
          sigma <<- s / fData$m
        }
        else{
          sigma[k] <<- sk / nk
        }
      }

    }
  )
)
