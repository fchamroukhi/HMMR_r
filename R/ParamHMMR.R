#' A Reference Class which contains parameters of a HMMR model.
#'
#' ParamHMMR contains all the parameters of a HMMR model.
#'
#' @field fData [FData][FData] object representing the sample.
#' @field K The number of regimes (mixture components).
#' @field p The order of the polynomial regression.
#' @field variance_type Character indicating if the model is homoskedastic
#' (`variance_type` = "homoskedastic") or heteroskedastic
#' (`variance_type` = "heteroskedastic").
#' @field prior The prior probabilities of the Markov chain.
#' @field trans_mat The transition matrix of the Markov chain.
#' @field beta Parameters of the polynomial regressions.
#' \eqn{\beta = (\beta_{1},\dots,\beta_{K})} is a matrix of dimension
#' \eqn{(p + 1, K)}, with \emph{p} the order of the polynomial regression.
#' @field sigma2 The variances for the \emph{K} regimes. If HMMR model is
#' homoskedastic (`variance_type` = "homoskedastic") then `sigma2` is a
#' matrix of size \eqn{(1, 1)}, else if HMMR model is heteroskedastic
#' (`variance_type` = "heteroskedastic") then `sigma2` is a matrix of size
#' \eqn{(K, 1)}.
#' @field nu The degree of freedom of the HMMR model.
#' @field phi A designed matrix for the polynomial regressions.
#' @seealso [FData]
#' @export
ParamHMMR <- setRefClass(
  "ParamHMMR",
  fields = list(
    fData = "FData",
    phi = "matrix",

    K = "numeric", # Number of regimes
    p = "numeric", # Dimension of beta (order of polynomial regression)
    variance_type = "character",
    nu = "numeric", # Degree of freedom

    prior = "matrix",
    trans_mat = "matrix",
    beta = "matrix",
    sigma2 = "matrix",
    mask = "matrix"
  ),
  methods = list(

    initialize = function(fData = FData(numeric(1), matrix(1)), K = 2, p = 2, variance_type = "heteroskedastic") {

      fData <<- fData

      phi <<- designmatrix(x = fData$X, p = p)$XBeta

      K <<- K
      p <<- p
      variance_type <<- variance_type

      if (variance_type == "homoskedastic") {
        nu <<- K - 1 + K * (K - 1) + K * (p + 1) + 1
      }
      else{
        nu <<- K - 1 + K * (K - 1) + K * (p + 1) + K
      }

      prior <<- matrix(NA, ncol = K)
      trans_mat <<- matrix(NA, K, K)
      beta <<- matrix(NA, p + 1, K)
      if (variance_type == "homoskedastic") {
        sigma2 <<- matrix(NA)
      }
      else{
        sigma2 <<- matrix(NA, K)
      }
      mask <<- matrix(NA, K, K)

    },

    initParam = function(try_algo = 1) {
      "Method to initialize parameters \\code{prior}, \\code{trans_mat},
      \\code{beta} and \\code{sigma2}.

      If try_algo = 1 then \\code{beta} and \\code{sigma2} are
      initialized by segmenting uniformly into \\code{K} contiguous segments
      the response Y. Otherwise, \\code{beta} and \\code{sigma2} are
      initialized by segmenting randomly the response Y into \\code{K} segments."

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
      if (try_algo == 1) { # Uniform segmentation into K contiguous segments, and then a regression

        zi <- round(fData$m / K) - 1

        s <- 0 # If homoskedastic
        for (k in 1:K) {
          yk <- fData$Y[((k - 1) * zi + 1):(k * zi)]
          Xk <- as.matrix(phi[((k - 1) * zi + 1):(k * zi), ])

          beta[, k] <<- solve(t(Xk) %*% Xk + (10 ^ -4) * diag(p + 1)) %*% t(Xk) %*% yk

          muk <- Xk %*% beta[, k]
          sk <- t(yk - muk) %*% (yk - muk)
          if (variance_type == "homoskedastic") {
            s <- (s + sk)
            sigma2 <<- s / fData$m
          }
          else {
            sigma2[k] <<- sk / length(yk)
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

          if (variance_type == "homoskedastic") {
            s <- s + sk
            sigma2 <<- s / fData$m

          }
          else{
            sigma2[k] <<- sk / length(yk)
          }
        }
      }

    },

    MStep = function(statHMMR) {
      "Method used in the EM algorithm to learn the parameters of the HMMR model
      based on statistics provided by \\code{statHMMR}."
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
        if (variance_type == "homoskedastic") {
          s <- (s + sk)
          sigma2 <<- s / fData$m
        }
        else{
          sigma2[k] <<- sk / nk
        }
      }

    }
  )
)
