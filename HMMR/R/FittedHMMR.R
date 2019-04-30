FittedHMMR <- setRefClass(
  "FittedHMMR",
  fields = list(
    modelHMMR = "ModelHMMR",
    paramHMMR = "ParamHMMR",
    statHMMR = "StatHMMR"
  ),
  methods = list(
    plot = function() {

      yaxislim <- c(mean(modelHMMR$Y) - 2 * sd(modelHMMR$Y), mean(modelHMMR$Y) + 2 * sd(modelHMMR$Y))

      # Predicted time series and predicted regime probabilities
      par(mfrow = c(2, 1), mai = c(0.6, 1, 0.5, 0.5), mgp = c(2, 1, 0))
      plot.default(modelHMMR$Y, type = "l", ylim = yaxislim, xlab = "x", ylab = "y")
      lines(statHMMR$predicted, type = "l", col = "red", lwd = 1.5)
      title(main = "Original and predicted HMMR time series")

      # Prediction probabilities of the hidden process (segmentation)
      colorsvec <- rainbow(modelHMMR$K)
      plot.default(statHMMR$predict_prob[, 1], type = "l", ylab = expression('P(Z'[t] == k ~ '|' ~ list(y[1],..., y[t - 1]) ~ ')'), col = colorsvec[1], lwd = 1.5, main = "Prediction probabilities")
      for (k in 2:modelHMMR$K) {
        lines(statHMMR$predict_prob[, k], col = colorsvec[k], lwd = 1.5) # Post Probs: Pr(Z_{t}=k|y_1,\ldots,y_{t-1})
      }

      # Filtered time series and filtering regime probabilities
      par(mfrow = c(2, 1), mai = c(0.6, 1, 0.5, 0.5), mgp = c(2, 1, 0))
      plot.default(modelHMMR$Y, type = "l", ylim = yaxislim, xlab = "x", ylab = "y")
      title(main = "Original and filtered HMMR time series")
      lines(statHMMR$filtered, col = "red", lwd = 1.5)

      # Filtering probabilities of the hidden process (segmentation)
      plot.default(statHMMR$filter_prob[, 1], type = "l", xlab = "t", ylab = expression('P(Z'[t] == k ~ '|' ~ list(y[1],..., y[t]) ~ ')'), col = colorsvec[1], lwd = 1.5)
      title(main = "Filtering probabilities")
      for (k in 2:modelHMMR$K) {
        lines(statHMMR$filter_prob[, k], col = colorsvec[k], lwd = 1.5) # Post Probs: Pr(Z_{t}=k|y_1,\ldots,y_t)
      }

      # Data, regressors, and segmentation
      par(mfrow = c(2, 1), mai = c(0.6, 1, 0.5, 0.5), mgp = c(2, 1, 0))
      plot.default(modelHMMR$Y, type = "l", ylim = yaxislim, xlab = "x", ylab = "y")
      title(main = "Time series, HMMR regimes, and smoothing probabilites")
      for (k in 1:modelHMMR$K) {
        model_k <- statHMMR$regressors[, k]
        #prob_model_k = HMMR$param$piik[,k]

        index <- statHMMR$klas == k
        active_model_k <- model_k[index] # prob_model_k >= prob);
        active_period_model_k <- seq(1:modelHMMR$m)[index] # prob_model_k >= prob);

        if (length(active_model_k) != 0) {
          lines(model_k, col = colorsvec[k], lty = "dotted", lwd = 1.5)
          lines(active_period_model_k, active_model_k, col = colorsvec[k], lwd = 1.5)
        }
      }

      # Probablities of the hidden process (segmentation)
      plot.default(statHMMR$tau_tk[, 1], type = "l", xlab = "x", ylab = expression('P(Z'[t] == k ~ '|' ~ list(y[1],..., y[n]) ~ ')'), col = colorsvec[1], lwd = 1.5)
      title(main = "Smoothing probabilities")
      if (modelHMMR$K > 1) {
        for (k in 2:modelHMMR$K) {
          lines(statHMMR$tau_tk[, k], col = colorsvec[k], lwd = 1.5) # Post Probs: Pr(Z_{t}=k|y_1,\ldots,y_n)
        }
      }

      # Data, regression model, and segmentation
      par(mfrow = c(2, 1), mai = c(0.6, 1, 0.5, 0.5), mgp = c(2, 1, 0))
      plot.default(modelHMMR$Y, type = "l", ylim = yaxislim, xlab = "x", ylab = "y")
      title(main = "Original and smoothed HMMR time series, and segmentation")
      lines(statHMMR$smoothed, col = "red", lwd = 1.5)

      # Transition time points
      tk <- which(diff(statHMMR$klas) != 0)
      for (i in 1:length(tk)) {
        abline(v = tk[i], col = "red", lty = "dotted", lwd = 2)
      }

      # Probablities of the hidden process (segmentation)
      plot.default(statHMMR$klas, type = "l", xlab = "x", ylab = "Estimated class labels", col = "red", lwd = 1.5)
    }
  )
)

FittedHMMR <- function(modelHMMR, paramHMMR, statHMMR) {
  new("FittedHMMR", modelHMMR = modelHMMR, paramHMMR = paramHMMR, statHMMR = statHMMR)
}
