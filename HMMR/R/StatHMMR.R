source("R/forwards_backwards.R")
source("R/hmm_process.R")

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
    cputime_total = "list", # cputime_total: for all the EM runs
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
    #           X: [nx(p+1)] regression design matrix
    #           nu: model complexity
    #           parameter_vector
  ),
  methods = list(
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
    # #######
    # # compute loglikelihood
    # #######
    computeLikelihood = function(paramHMMR) {
      fb <- forwards_backwards(paramHMMR$prior, paramHMMR$trans_mat, f_tk)
      loglik <<- fb$loglik

    },
    #######
    # compute the final solution stats
    #######
    computeStats = function(modelHMMR, paramHMMR, phi, cputime_total) {
      cputime <<- mean(cputime_total)

      # ## sate sequence prob p(z_1,...,z_n;\pi,A)
      state_probs <<- hmm_process(paramHMMR$prior, paramHMMR$trans_mat, modelHMMR$m)

      ### BIC, AIC, ICL
      BIC <<- loglik - modelHMMR$nu * log(modelHMMR$m) / 2
      AIC <<- loglik - modelHMMR$nu

      # # CL(theta) : Completed-data loglikelihood
      # sum_t_log_Pz_ftk = sum(hmmr.stats.Zik.*log(state_probs.*hmmr.stats.f_tk), 2);
      # comp_loglik = sum(sum_t_log_Pz_ftk(K:end));
      # hmmr.stats.comp_loglik = comp_loglik;
      # hmmr.stats.ICL = comp_loglik - (nu*log(m)/2);

      ## predicted, filtered, and smoothed time series
      regressors <<- round(phi %*% paramHMMR$beta, 4)

      # prediction probs   = Pr(z_t|y_1,...,y_{t-1})
      predict_prob[1, ] <<- paramHMMR$prior # t=1 p (z_1)
      predict_prob[2:modelHMMR$m, ] <<- round((alpha_tk[(1:(modelHMMR$m - 1)), ] %*% paramHMMR$trans_mat) / (apply(alpha_tk[(1:(modelHMMR$m - 1)), ], 1, sum) %*% matrix(1, 1, modelHMMR$K)), 5) #t =2,...,n

      # predicted observations
      predicted <<- matrix(apply(round(predict_prob * regressors, 5), 1, sum)) #pond par les probas de prediction

      # filtering probs  = Pr(z_t|y_1,...,y_t)
      filter_prob <<- round(alpha_tk / (apply(alpha_tk, 1, sum) %*% matrix(1, 1, modelHMMR$K)), 5) #normalize(alpha_tk,2);

      # filetered observations
      filtered <<- as.matrix(apply(round(filter_prob * regressors, 5), 1, sum)) #pond par les probas de filtrage

      ### smoothed observations
      smoothed_regressors <<- tau_tk * regressors
      smoothed <<- as.matrix(apply(smoothed_regressors, 1, sum))

    },
    #######
    # EStep
    #######
    EStep = function(modelHMMR, paramHMMR, phi) {
      muk <- matrix(0, modelHMMR$m, modelHMMR$K)

      # observation likelihoods
      for (k in 1:modelHMMR$K) {
        mk <- phi %*% paramHMMR$beta[, k]
        muk[, k] <- mk
        # the regressors means
        if (modelHMMR$variance_type == variance_types$homoskedastic) {
          sk <- paramHMMR$sigma
        }
        else{
          sk <- paramHMMR$sigma[k]
        }
        z <- ((modelHMMR$Y - mk) ^ 2) / sk
        log_f_tk[, k] <<- -0.5 * matrix(1, modelHMMR$m, 1) %*% (log(2 * pi) + log(sk)) - 0.5 * z#log(gaussienne)

      }

      log_f_tk <<- pmin(log_f_tk, log(.Machine$double.xmax))
      log_f_tk <<- pmax(log_f_tk, log(.Machine$double.xmin))

      f_tk <<- exp(log_f_tk)

      fb <- forwards_backwards(paramHMMR$prior, paramHMMR$trans_mat, f_tk)
      tau_tk <<- fb$tau_tk
      xi_tkl <<- fb$xi_tkl
      alpha_tk <<- fb$alpha_tk
      beta_tk <<- fb$beta_tk
      loglik <<- fb$loglik

    }
  )
)


StatHMMR <- function(modelHMMR) {
  tau_tk <- matrix(NA, modelHMMR$m, modelHMMR$K) # tau_tk: smoothing probs: [nxK], tau_tk(t,k) = Pr(z_i=k | y1...yn)
  alpha_tk <- matrix(NA, modelHMMR$m, ncol = modelHMMR$K) # alpha_tk: [nxK], forwards probs: Pr(y1...yt,zt=k)
  beta_tk <- matrix(NA, modelHMMR$m, modelHMMR$K) # beta_tk: [nxK], backwards probs: Pr(yt+1...yn|zt=k)
  xi_tkl <- array(NA, c(modelHMMR$m - 1, modelHMMR$K, modelHMMR$K)) # xi_tkl: [(n-1)xKxK], joint post probs : xi_tk\elll(t,k,\ell)  = Pr(z_t=k, z_{t-1}=\ell | Y) t =2,..,n
  f_tk <- matrix(NA, modelHMMR$m, modelHMMR$K) # f_tk: [nxK] f(yt|zt=k)
  log_f_tk <- matrix(NA, modelHMMR$m, modelHMMR$K) # log_f_tk: [nxK] log(f(yt|zt=k))
  loglik <- -Inf # loglik: log-likelihood at convergence
  stored_loglik <- list() # stored_loglik: stored log-likelihood values during EM
  cputime <- Inf # cputime: for the best run
  cputime_total <- list() # cputime_total: for all the EM runs
  klas <- matrix(NA, modelHMMR$m, 1) # klas: [nx1 double]
  z_ik <- matrix(NA, modelHMMR$m, modelHMMR$K) # z_ik: [nxK]
  state_probs <- matrix(NA, modelHMMR$m, modelHMMR$K) # state_probs: [nxK]
  BIC <- -Inf # BIC
  AIC <- -Inf # AIC
  regressors <- matrix(NA, modelHMMR$m, modelHMMR$K) # regressors: [nxK]
  predict_prob <- matrix(NA, modelHMMR$m, modelHMMR$K) # predict_prob: [nxK]: Pr(zt=k|y1...y_{t-1})
  predicted <- matrix(NA, modelHMMR$m, 1) # predicted: [nx1]
  filter_prob <- matrix(NA, modelHMMR$m, modelHMMR$K) # filter_prob: [nxK]: Pr(zt=k|y1...y_t)
  filtered <- matrix(NA, modelHMMR$m, 1) # filtered: [nx1]
  smoothed_regressors <- matrix(NA, modelHMMR$m, modelHMMR$K) # smoothed_regressors: [nxK]
  smoothed <- matrix(NA, modelHMMR$m, 1) # smoothed: [nx1]

  new(
    "StatHMMR",
    tau_tk = tau_tk,
    alpha_tk = alpha_tk,
    beta_tk = beta_tk,
    xi_tkl = xi_tkl,
    f_tk = f_tk,
    log_f_tk = log_f_tk,
    loglik = loglik,
    stored_loglik = stored_loglik,
    cputime = cputime,
    cputime_total = cputime_total,
    klas = klas,
    z_ik = z_ik,
    state_probs = state_probs,
    BIC = BIC,
    AIC = AIC,
    regressors = regressors,
    predict_prob = predict_prob,
    predicted = predicted,
    filter_prob = filter_prob,
    filtered = filtered,
    smoothed_regressors = smoothed_regressors,
    smoothed = smoothed
  )
}
