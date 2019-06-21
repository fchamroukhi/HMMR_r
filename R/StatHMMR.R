#' A Reference Class which contains statistics of a HMMR model.
#'
#' StatHMMR contains all the parameters of a [HMMR][ParamHMMR] model.
#'
#' @field tau_tk Matrix of size \eqn{(m, K)} giving the posterior probability that
#' the observation \eqn{Y_{t}}{Y_{t}} originates from the \eqn{k}-th
#' regression model \eqn{P(z\_t = k | Y_{1},\dots,Y_{m})}{P(z_t = k | Y_{1},\dots,Y_{m})}.
#' @field alpha_tk Matrix of size \eqn{(m, K)} giving the forwards
#' probabilities: \eqn{P(Y_{1},\dots,Y_{t}, z\_t = k)}{P(Y_{1},\dots,Y_{t}, z_t = k)}.
#' @field beta_tk Matrix of size \eqn{(m, K)}, giving the backwards
#' probabilities: \eqn{P(Y_{t+1},\dots,Y_{m} | z\_t = k)}{P(Y_{t+1},\dots,Y_{m} | z_t = k)}.
#' @field xi_tkl Array of size \eqn{(m - 1, K, K)} giving the joint post
#' probabilities: \eqn{xi_tk[t, k, l] = P(z\_t = k, z_{t-1} = l | Y)}{xi_tk[t, k, l] = P(z_t = k, z_{t-1} = l | Y)}
#' for \eqn{t = 2,\dots,m}.
#' @field f_tk Matrix of size \eqn{(m, K)} giving the cumulative distribution
#' function \eqn{f(yt | z\_t = k)}{f(yt | z_t = k)}.
#' @field log_f_tk Matrix of size \eqn{(m, K)} giving the logarithm of the
#' cumulative distribution `f_tk`.
#' @field loglik Numeric. Log-likelihood of the HMMR model.
#' @field stored_loglik List. Stored values of the log-likelihood at each
#' iteration of the EM algorithm.
#' @field cpu_time Numeric. Average executon time of the EM algorithm. for the best run.
#' @field klas Column matrix of the labels issued from `z_ik`. Its elements are
#' \eqn{klas(i) = k}, \eqn{k = 1,\dots,K}.
#' @field z_ik Hard segmentation logical matrix of dimension \eqn{(m, K)}
#' obtained by the Maximum a posteriori (MAP) rule:
#' \eqn{z_{ik} = 1 \ \textrm{if} \ z_{ik} = \textrm{arg} \ \textrm{max}_{k} \
#' P(z_i = k | Y) = tau_ik;\ 0 \ \textrm{otherwise}}{z_ik = 1 if z_ik =
#' arg max_k P(z_i = k | Y) = tau_ik; 0 otherwise}, \eqn{k = 1,\dots,K}.
#' @field state_probs Matrix of size \eqn{(m, K)} giving the distribution of
#' the Markov chain \eqn{P(z_{1},\dots,z_{m};\pi,A)}{P(z_{1},\dots,z_{m};\pi,A)}
#' with \eqn{\pi} the prior probabilities (field `prior` of the class
#' [ParamHMMR][ParamHMMR]) and \eqn{A} the transition matrix (field `trans_mat`
#' of the class [ParamHMMR][ParamHMMR]) of the Markov chain.
#' @field BIC Numeric. Value of the BIC (Bayesian Information Criterion)
#' criterion. The formula is \eqn{BIC = loglik - nu \times
#' \textrm{log}(m) / 2}{BIC = loglik - nu x log(m) / 2} with `nu` the
#' degree of freedom of the HMMR model.
#' @field AIC Numeric. Value of the AIC (Akaike Information Criterion)
#' criterion. The formula is \eqn{AIC = loglik - nu}.
#' @field regressors Matrix of size \eqn{(m, K)} giving the values of
#' \eqn{\beta_{k} \times X_{i}}{\beta_{k} x X_{i}}, \eqn{i = 1,\dots,m}.
#' @field predict_prob Matrix of size \eqn{(m, K)} giving the prediction
#' probabilities: \eqn{P(z\_t = k | y_{1},\dots,y_{t-1})}{P(z_t = k | y_{1},\dots,y_{t-1})}.
#' @field predicted Row matrix of size \eqn{(m, 1)} giving the predicted
#' observations \eqn{\sum_{k} P(z\_t = k | y_{1},\dots,y_{t-1}) \times \beta_{k} \times X_{i}}{\sum_{k} P(z_t = k | y_{1},\dots,y_{t-1}) x \betak x X_{i}}.
#' @field filter_prob Matrix of size \eqn{(m, K)} giving the filtering
#' probabilities \eqn{Pr(z\_t = k | y_{1},\dots,y_{t})}{Pr(z_t = k | y_{1},\dots,y_{t})}.
#' @field filtered Row matrix of size \eqn{(m, 1)} giving the filetered
#' observations \eqn{\sum_{k} P(z\_t = k | y_{1},\dots,y_{t}) \times \beta_{k} \times X_{i}}{\sum_{k} P(z_t = k | y_{1},\dots,y_{t}) x \betak x X_{i}}.
#' @field smoothed_regressors Matrix of size \eqn{(m, K)} giving the smoothed observations:
#' \eqn{P(z\_i = k | Y_{1},\dots,Y_{n}) \times \beta_{k} \times X_{i}}{P(z_i = k | Y_{1},\dots,Y_{n}) x \betak x X_{i}}.
#' @field smoothed Row matrix of size \eqn{(m, 1)} giving:
#' \eqn{\sum_{k} P(z\_i = k | Y_{1},\dots,Y_{n}) \times \beta_{k} \times X_{i}}{\sum_{k}  P(z_i = k | Y_{1},\dots,Y_{n}) x \betak x X_{i}}
#' @seealso [ParamHMMR]
#' @export
StatHMMR <- setRefClass(
  "StatHMMR",
  fields = list(
    tau_tk = "matrix", # tau_tk: smoothing probs: [nxK], tau_tk(t,k) = Pr(z_i=k | y1...yn)
    alpha_tk = "matrix", # alpha_tk: [nxK], forwards probs: Pr(y1...yt,zt=k)
    beta_tk = "matrix", # beta_tk: [nxK], backwards probs: Pr(yt+1...yn|zt=k)
    xi_tkl = "array", # xi_tkl: [(n-1)xKxK], joint post probs : xi_tk\elll(t,k,\ell)  = Pr(z_t=k, z_{t-1}=\ell | Y) t =2,..,n
    f_tk = "matrix", # f_tk: [nxK] f(yt|zt=k)
    log_f_tk = "matrix", # log_f_tk: [nxK] log(f(yt|zt=k))
    loglik = "numeric", # loglik: log-likelihood at convergence
    stored_loglik = "list", # stored_loglik: stored log-likelihood values during EM
    cputime = "numeric", # cputime: for the best run
    klas = "matrix", # klas: [nx1 double]
    z_ik = "matrix", # z_ik: [nxK]
    state_probs = "matrix", # state_probs: [nxK]
    BIC = "numeric", # BIC
    AIC = "numeric", # AIC
    regressors = "matrix", # regressors: [nxK]
    predict_prob = "matrix", # predict_prob: [nxK]: Pr(zt=k|y1...y_{t-1})
    predicted = "matrix", # predicted: [nx1]
    filter_prob = "matrix", # filter_prob: [nxK]: Pr(zt=k|y1...y_t)
    filtered = "matrix", # filtered: [nx1]
    smoothed_regressors = "matrix", # smoothed_regressors: [nxK]
    smoothed = "matrix" # smoothed: [nx1]
  ),
  methods = list(

    initialize = function(paramHMMR = ParamHMMR()) {

      tau_tk <<- matrix(NA, paramHMMR$m, paramHMMR$K) # tau_tk: smoothing probs: [nxK], tau_tk(t,k) = Pr(z_i=k | y1...yn)
      alpha_tk <<- matrix(NA, paramHMMR$m, ncol = paramHMMR$K) # alpha_tk: [nxK], forwards probs: Pr(y1...yt,zt=k)
      beta_tk <<- matrix(NA, paramHMMR$m, paramHMMR$K) # beta_tk: [nxK], backwards probs: Pr(yt+1...yn|zt=k)
      xi_tkl <<- array(NA, c(paramHMMR$m - 1, paramHMMR$K, paramHMMR$K)) # xi_tkl: [(n-1)xKxK], joint post probs : xi_tk\elll(t,k,\ell)  = Pr(z_t=k, z_{t-1}=\ell | Y) t =2,..,n
      f_tk <<- matrix(NA, paramHMMR$m, paramHMMR$K) # f_tk: [nxK] f(yt|zt=k)
      log_f_tk <<- matrix(NA, paramHMMR$m, paramHMMR$K) # log_f_tk: [nxK] log(f(yt|zt=k))
      loglik <<- -Inf # loglik: log-likelihood at convergence
      stored_loglik <<- list() # stored_loglik: stored log-likelihood values during EM
      cputime <<- Inf # cputime: for the best run
      klas <<- matrix(NA, paramHMMR$m, 1) # klas: [nx1 double]
      z_ik <<- matrix(NA, paramHMMR$m, paramHMMR$K) # z_ik: [nxK]
      state_probs <<- matrix(NA, paramHMMR$m, paramHMMR$K) # state_probs: [nxK]
      BIC <<- -Inf # BIC
      AIC <<- -Inf # AIC
      regressors <<- matrix(NA, paramHMMR$m, paramHMMR$K) # regressors: [nxK]
      predict_prob <<- matrix(NA, paramHMMR$m, paramHMMR$K) # predict_prob: [nxK]: Pr(zt=k|y1...y_{t-1})
      predicted <<- matrix(NA, paramHMMR$m, 1) # predicted: [nx1]
      filter_prob <<- matrix(NA, paramHMMR$m, paramHMMR$K) # filter_prob: [nxK]: Pr(zt=k|y1...y_t)
      filtered <<- matrix(NA, paramHMMR$m, 1) # filtered: [nx1]
      smoothed_regressors <<- matrix(NA, paramHMMR$m, paramHMMR$K) # smoothed_regressors: [nxK]
      smoothed <<- matrix(NA, paramHMMR$m, 1) # smoothed: [nx1]

    },

    MAP = function() {
      N <- nrow(tau_tk)
      K <- ncol(tau_tk)
      ikmax <- max.col(tau_tk)
      ikmax <- matrix(ikmax, ncol = 1)
      z_ik <<- ikmax %*% ones(1, K) == ones(N, 1) %*% (1:K) # partition_MAP
      klas <<- ones(N, 1)
      for (k in 1:K) {
        klas[z_ik[, k] == 1] <<- k
      }
    },

    computeLikelihood = function(paramHMMR) {
      fb <- forwardsBackwards(paramHMMR$prior, paramHMMR$trans_mat, t(f_tk))
      loglik <<- fb$loglik

    },

    computeStats = function(paramHMMR, cputime_total) {
      cputime <<- mean(cputime_total)

      # State sequence prob p(z_1,...,z_n;\pi,A)
      state_probs <<- hmmProcess(paramHMMR$prior, paramHMMR$trans_mat, paramHMMR$m)

      # BIC, AIC, ICL
      BIC <<- loglik - paramHMMR$nu * log(paramHMMR$m) / 2
      AIC <<- loglik - paramHMMR$nu

      # # CL(theta) : Completed-data loglikelihood
      # sum_t_log_Pz_ftk = sum(hmmr.stats.Zik.*log(state_probs.*hmmr.stats.f_tk), 2);
      # comp_loglik = sum(sum_t_log_Pz_ftk(K:end));
      # hmmr.stats.comp_loglik = comp_loglik;
      # hmmr.stats.ICL = comp_loglik - (nu*log(m)/2);

      # Predicted, filtered, and smoothed time series
      regressors <<- paramHMMR$phi %*% paramHMMR$beta

      # Prediction probabilities = Pr(z_t|y_1,...,y_{t-1})
      predict_prob[1, ] <<- paramHMMR$prior # t=1 p (z_1)

      predict_prob[2:paramHMMR$m, ] <<- (alpha_tk[(1:(paramHMMR$m - 1)), ] %*% paramHMMR$trans_mat) / (apply(as.matrix(alpha_tk[(1:(paramHMMR$m - 1)), ]), 1, sum) %*% matrix(1, 1, paramHMMR$K)) # t = 2,...,n

      # Predicted observations
      predicted <<- matrix(apply(predict_prob * regressors, 1, sum)) # Weighted by prediction probabilities

      # Filtering probabilities = Pr(z_t|y_1,...,y_t)
      filter_prob <<- alpha_tk / (apply(alpha_tk, 1, sum) %*% matrix(1, 1, paramHMMR$K)) # Normalize(alpha_tk,2);

      # Filetered observations
      filtered <<- as.matrix(apply(filter_prob * regressors, 1, sum)) # Weighted by filtering probabilities

      # Smoothed observations
      smoothed_regressors <<- tau_tk * regressors
      smoothed <<- as.matrix(apply(smoothed_regressors, 1, sum))

    },

    EStep = function(paramHMMR) {
      muk <- matrix(0, paramHMMR$m, paramHMMR$K)

      # Observation likelihoods
      for (k in 1:paramHMMR$K) {
        mk <- paramHMMR$phi %*% paramHMMR$beta[, k]
        muk[, k] <- mk
        # The regressors means
        if (paramHMMR$variance_type == "homoskedastic") {
          sk <- paramHMMR$sigma2[1]
        }
        else{
          sk <- paramHMMR$sigma2[k]
        }
        z <- ((paramHMMR$Y - mk) ^ 2) / sk
        log_f_tk[, k] <<- -0.5 * matrix(1, paramHMMR$m, 1) %*% (log(2 * pi) + log(sk)) - 0.5 * z # log(gaussienne)

      }

      log_f_tk <<- pmin(log_f_tk, log(.Machine$double.xmax))
      log_f_tk <<- pmax(log_f_tk, log(.Machine$double.xmin))

      f_tk <<- exp(log_f_tk)

      fb <- forwardsBackwards(paramHMMR$prior, paramHMMR$trans_mat, t(f_tk))
      tau_tk <<- t(fb$tau_tk)
      xi_tkl <<- fb$xi_tkl
      alpha_tk <<- t(fb$alpha_tk)
      beta_tk <<- t(fb$beta_tk)
      loglik <<- fb$loglik

    }
  )
)
