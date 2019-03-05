show_HMMR_results <- function(x,y,HMMR,yaxislim){
#

if (nrow(x)!=1){
x = t(x); #y=y'
}
colors = c("red","blue","green","pink","cadetblue2","orange","blue4","chartreuse4","brown2","cadetblue4")

if ((nargs()<4)||length(yaxislim)==0){
  yaxislim = c(mean(y)-2*sd(y), mean(y)+2*sd(y))
}
K = ncol(hmmr$stats$tau_tk)

## predicted time series and predicted regime probabilities
x11()
par(mfrow=c(2,1))
plot(x,y,main="Original and predicted HMMR time series",ylab = 'y',xlab="",ylim=yaxislim,type="l")#black
lines(x,hmmr$stats$predicted,col="red",type="l",lwd=2)

## prediction probabilities of the hidden process (segmentation)
plot(x,hmmr$stats$predict_prob[,1],col=colors[1],type="l",lwd=1.5,main="prediction probabilities",xlab="t",ylab="Prob")
for (k in 2:K){
  lines(x,hmmr$stats$predict_prob[,k],col=colors[k]);
}

## filtered time series and filtering regime probabilities
x11()
par(mfrow=c(2,1))
plot(x,y,main="Original and filtered HMMR time series",type="l",ylab="y",xlab="",ylim=yaxislim)#black
lines(x,hmmr$stats$filtered,col="red",lwd=2)

## filtering probabilities of the hidden process (segmentation)
plot(x,hmmr$stats$filter_prob[,1],col=colors[1],type="l",lwd=1.5,main="filtering probabilities",xlab="t",ylab="Prob")
for (k in 2:K){
  lines(x,hmmr$stats$filter_prob[,k],col=colors[k],type="l",lwd=1.5) #Post Probs: Pr(Z_{t}=k|y_1,\ldots,y_t)
}

## data, regressors, and segmentation
x11()
par(mfrow=c(2,1))
plot(x,y,main="Time series, HMMR regimes, and smoothing probabilites",type="l",ylab="y",xlab="",ylim=yaxislim)
K = ncol(hmmr$stats$tau_tk)
for (k in 1:K){
  model_k = hmmr$stats$regressors[,k]
  #prob_model_k = HMMR$param$piik[,k]

  active_model_k = model_k[hmmr$stats$klas==k]#prob_model_k >= prob);
  active_period_model_k = x[hmmr$stats$klas==k]#prob_model_k >= prob);

  inactive_model_k = model_k[hmmr$stats$klas != k]#prob_model_k >= prob);
  inactive_period_model_k = x[hmmr$stats$klas != k]#prob_model_k >= prob);

  if (length(active_model_k)!=0){
    lines(inactive_period_model_k,inactive_model_k,col=colors[k],type="l",lwd=0.001)
    lines(active_period_model_k, active_model_k, col=colors[k],type="l",lwd=5)
  }
}

# # Probablities of the hidden process (segmentation)
plot(x,hmmr$stats$tau_tk[,1],main="smoothing probabilities",xlab="t",ylab="Prob",col=colors[1],type="l",lwd=1.5)
for (k in 2:K){
  lines(x,hmmr$stats$tau_tk[,k],col=colors[k],type="l",lwd=1.5)
}
#Post Probs: Pr(Z_{t}=k|y_1,\ldots,y_n)
 
## data, regression model, and segmentation
x11()
par(mfrow=c(2,1))
plot(x,y,main="Original and smoothed HMMR time series, and segmentation",type="l",ylab="y",xlab="",ylim=yaxislim)#black
lines(x,hmmr$stats$smoothed,col="red",lwd=2)

# transition time points
tk = which(diff(hmmr$stats$klas)!=0)
lines(rbind(x[tk],x[tk]),t(matrix(c(1),length(tk),1)%*%c(min(y)-2*sd(y), max(y)+2*sd(y))),type='h',lty=2,lwd=2)

##Probablities of the hidden process (segmentation)
plot(x,hmmr$stats$klas,lwd=1.5, xlab="t",ylab="Estimated class labels")

#### model log-likelihood during EM
# plot(hmmr$stats$stored_loglik,type='l',lwd=1.5, xlab="EM iteration number",ylab="log-likelihood")

}
