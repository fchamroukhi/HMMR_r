#' @export
ModelHMMR <- setRefClass(
  "ModelHMMR",
  fields = list(
    paramHMMR = "ParamHMMR",
    statHMMR = "StatHMMR"
  ),
  methods = list(
    plot = function() {

      oldpar <- par()[c("mfrow", "mai", "mgp")]
      on.exit(par(oldpar), add = TRUE)

      yaxislim <- c(mean(paramHMMR$fData$Y) - 2 * sd(paramHMMR$fData$Y), mean(paramHMMR$fData$Y) + 2 * sd(paramHMMR$fData$Y))

      # Predicted time series and predicted regime probabilities
      par(mfrow = c(2, 1), mai = c(0.6, 1, 0.5, 0.5), mgp = c(2, 1, 0))
      plot.default(paramHMMR$fData$X, paramHMMR$fData$Y, type = "l", ylim = yaxislim, xlab = "x", ylab = "y")
      lines(paramHMMR$fData$X, statHMMR$predicted, type = "l", col = "red", lwd = 1.5)
      title(main = "Original and predicted HMMR time series")

      # Prediction probabilities of the hidden process (segmentation)
      colorsvec <- rainbow(paramHMMR$K)
      plot.default(paramHMMR$fData$X, statHMMR$predict_prob[, 1], type = "l", xlab = "x", ylab = expression('P(Z'[t] == k ~ '|' ~ list(y[1],..., y[t - 1]) ~ ')'), col = colorsvec[1], lwd = 1.5, main = "Prediction probabilities", ylim = c(0, 1))
      for (k in 2:paramHMMR$K) {
        lines(paramHMMR$fData$X, statHMMR$predict_prob[, k], col = colorsvec[k], lwd = 1.5) # Post Probs: Pr(Z_{t}=k|y_1,\ldots,y_{t-1})
      }

      # Filtered time series and filtering regime probabilities
      par(mfrow = c(2, 1), mai = c(0.6, 1, 0.5, 0.5), mgp = c(2, 1, 0))
      plot.default(paramHMMR$fData$X, paramHMMR$fData$Y, type = "l", ylim = yaxislim, xlab = "x", ylab = "y")
      title(main = "Original and filtered HMMR time series")
      lines(paramHMMR$fData$X, statHMMR$filtered, col = "red", lwd = 1.5)

      # Filtering probabilities of the hidden process (segmentation)
      plot.default(paramHMMR$fData$X, statHMMR$filter_prob[, 1], type = "l", xlab = "x", ylab = expression('P(Z'[t] == k ~ '|' ~ list(y[1],..., y[t]) ~ ')'), col = colorsvec[1], lwd = 1.5, ylim = c(0, 1))
      title(main = "Filtering probabilities")
      for (k in 2:paramHMMR$K) {
        lines(paramHMMR$fData$X, statHMMR$filter_prob[, k], col = colorsvec[k], lwd = 1.5) # Post Probs: Pr(Z_{t}=k|y_1,\ldots,y_t)
      }

      # Data, regressors, and segmentation
      par(mfrow = c(2, 1), mai = c(0.6, 1, 0.5, 0.5), mgp = c(2, 1, 0))
      plot.default(paramHMMR$fData$X, paramHMMR$fData$Y, type = "l", ylim = yaxislim, xlab = "x", ylab = "y")
      title(main = "Time series, HMMR regimes, and smoothing probabilites")
      for (k in 1:paramHMMR$K) {
        model_k <- statHMMR$regressors[, k]
        #prob_model_k = HMMR$param$piik[,k]

        index <- statHMMR$klas == k
        active_model_k <- model_k[index] # prob_model_k >= prob);
        active_period_model_k <- paramHMMR$fData$X[index] # prob_model_k >= prob);

        if (length(active_model_k) != 0) {
          lines(paramHMMR$fData$X, model_k, col = colorsvec[k], lty = "dotted", lwd = 1.5)
          lines(active_period_model_k, active_model_k, col = colorsvec[k], lwd = 1.5)
        }
      }

      # Probablities of the hidden process (segmentation)
      plot.default(paramHMMR$fData$X, statHMMR$tau_tk[, 1], type = "l", xlab = "x", ylab = expression('P(Z'[t] == k ~ '|' ~ list(y[1],..., y[n]) ~ ')'), col = colorsvec[1], lwd = 1.5, ylim = c(0, 1))
      title(main = "Smoothing probabilities")
      if (paramHMMR$K > 1) {
        for (k in 2:paramHMMR$K) {
          lines(paramHMMR$fData$X, statHMMR$tau_tk[, k], col = colorsvec[k], lwd = 1.5) # Post Probs: Pr(Z_{t}=k|y_1,\ldots,y_n)
        }
      }

      # Data, regression model, and segmentation
      par(mfrow = c(2, 1), mai = c(0.6, 1, 0.5, 0.5), mgp = c(2, 1, 0))
      plot.default(paramHMMR$fData$X, paramHMMR$fData$Y, type = "l", ylim = yaxislim, xlab = "x", ylab = "y")
      title(main = "Original and smoothed HMMR time series, and segmentation")
      lines(paramHMMR$fData$X, statHMMR$smoothed, col = "red", lwd = 1.5)

      # Transition time points
      tk <- which(diff(statHMMR$klas) != 0)
      for (i in 1:length(tk)) {
        abline(v = tk[i], col = "red", lty = "dotted", lwd = 2)
      }

      # Probablities of the hidden process (segmentation)
      plot.default(paramHMMR$fData$X, statHMMR$klas, type = "l", xlab = "x", ylab = "Estimated class labels", col = "red", lwd = 1.5)
      axis(side = 2, at = 1:paramHMMR$K)
    },

    summary = function() {

      digits = getOption("digits")

      title <- paste("Fitted HMMR model")
      txt <- paste(rep("-", min(nchar(title) + 4, getOption("width"))), collapse = "")

      # Title
      cat(txt)
      cat("\n")
      cat(title)
      cat("\n")
      cat(txt)

      cat("\n")
      cat("\n")
      cat(paste0("HMMR model with K = ", paramHMMR$K, ifelse(paramHMMR$K > 1, " components", " component"), ":"))
      cat("\n")
      cat("\n")

      tab <- data.frame("log-likelihood" = statHMMR$loglik, "nu" = paramHMMR$nu, "AIC" = statHMMR$AIC,
                        "BIC" = statHMMR$BIC, row.names = "", check.names = FALSE)
      print(tab, digits = digits)

      cat("\nClustering table (Number of observations in each regimes):\n")
      print(table(statHMMR$klas))

      cat("\nRegression coefficients:\n\n")
      if (paramHMMR$p > 0) {
        row.names = c("1", sapply(1:paramHMMR$p, function(x) paste0("X^", x)))
      } else {
        row.names = "1"
      }

      betas <- data.frame(paramHMMR$beta, row.names = row.names)
      colnames(betas) <- sapply(1:paramHMMR$K, function(x) paste0("Beta(K = ", x, ")"))
      print(betas, digits = digits)

      cat(paste0(ifelse(paramHMMR$variance_type == "homoskedastic", "\n\n",
                        "\nVariances:\n\n")))
      sigma2 = data.frame(t(paramHMMR$sigma2), row.names = NULL)
      if (paramHMMR$variance_type == "homoskedastic") {
        colnames(sigma2) = "Sigma2"
        print(sigma2, digits = digits, row.names = FALSE)
      } else {
        colnames(sigma2) = sapply(1:paramHMMR$K, function(x) paste0("Sigma2(K = ", x, ")"))
        print(sigma2, digits = digits, row.names = FALSE)
      }

    }
  )
)
