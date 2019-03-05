learn_hmmr<- function(x, y, K, p,type_variance, total_EM_tries, max_iter_EM, threshold, verbose){
  
  # learn_hmmr learn a Regression model with a Hidden Markov Process (HMMR)
  # for modeling and segmentation of a time series with regime changes.
  # The learning is performed by the EM (Baum-Welch) algorithm.
  #
  #
  # Inputs :
  #
  #          (x,y) : a time series composed of m points : dim(y)=[m 1]
  #                * Each curve is observed during the interval [0,T], i.e x =[t_1,...,t_m]
  #
  #           K : Number of polynomial regression components (regimes)
  #          	p : degree of the polynomials
  #
  # Outputs :
  #
  #         hmmr: the estimated HMMR model. a structure composed of:
  #
  #         prior: [Kx1]: prior(k) = Pr(z_1=k), k=1...K
  #         trans_mat: [KxK], trans_mat(\ell,k) = Pr(z_t = k|z_{t-1}=\ell)
  #         reg_param: the paramters of the regressors:
  #                 betak: regression coefficients
  #                 sigma2k (or sigma2) : the variance(s)
  #         Stats:
  #           tau_tk: smoothing probs: [nxK], tau_tk(t,k) = Pr(z_i=k | y1...yn)
  #           alpha_tk: [nxK], forwards probs: Pr(y1...yt,zt=k)
  #           beta_tk: [nxK], backwards probs: Pr(yt+1...yn|zt=k)
  #           xi_tkl: [(n-1)xKxK], joint post probs : xi_tk\elll(t,k,\ell)  = Pr(z_t=k, z_{t-1}=\ell | Y) t =2,..,n
  #           X: [nx(p+1)] regression design matrix
  #           nu: model complexity
  #           parameter_vector
  #           f_tk: [nxK] f(yt|zt=k)
  #           log_f_tk: [nxK] log(f(yt|zt=k))
  #           loglik: log-likelihood at convergence
  #           stored_loglik: stored log-likelihood values during EM
  #           cputime: for the best run
  #           cputime_total: for all the EM runs
  #           klas: [nx1 double]
  #           Zik: [nxK]
  #           state_probs: [nxK]
  #           BIC: -2.1416e+03
  #           AIC: -2.0355e+03
  #           regressors: [nxK]
  #           predict_prob: [nxK]: Pr(zt=k|y1...y_{t-1})
  #           predicted: [nx1]
  #           filter_prob: [nxK]: Pr(zt=k|y1...y_t)
  #           filtered: [nx1]
  #           smoothed_regressors: [nxK]
  #           smoothed: [nx1]
  #
  #
  #Faicel Chamroukhi, sept 2008 
  #
  ## Please cite the following papers for this code:
  #
  # 
  # @article{Chamroukhi-FDA-2018,
  # 	Journal = {Wiley Interdisciplinary Reviews: Data Mining and Knowledge Discovery},
  # 	Author = {Faicel Chamroukhi and Hien D. Nguyen},
  # 	Note = {DOI: 10.1002/widm.1298.},
  # 	Volume = {},
  # 	Title = {Model-Based Clustering and Classification of Functional Data},
  # 	Year = {2019},
  # 	Month = {to appear},
  # 	url =  {https://chamroukhi.com/papers/MBCC-FDA.pdf}
  # 	}
  # 
  # @InProceedings{Chamroukhi-IJCNN-2011,
  #   author = {F. Chamroukhi and A. Sam\'e  and P. Aknin and G. Govaert},
  #   title = {Model-based clustering with Hidden Markov Model regression for time series with regime changes},
  #   Booktitle = {Proceedings of the International Joint Conference on Neural Networks (IJCNN), IEEE},
  #   Pages = {2814--2821},
  #   Adress = {San Jose, California, USA},
  #   year = {2011},
  #   month = {Jul-Aug},
  #   url = {https://chamroukhi.com/papers/Chamroukhi-ijcnn-2011.pdf}
  # }
  # 
  # @INPROCEEDINGS{Chamroukhi-IJCNN-2009,
  #   AUTHOR =       {Chamroukhi, F. and Sam\'e,  A. and Govaert, G. and Aknin, P.},
  #   TITLE =        {A regression model with a hidden logistic process for feature extraction from time series},
  #   BOOKTITLE =    {International Joint Conference on Neural Networks (IJCNN)},
  #   YEAR =         {2009},
  #   month = {June},
  #   pages = {489--496},
  #   Address = {Atlanta, GA},
  #  url = {https://chamroukhi.com/papers/chamroukhi_ijcnn2009.pdf},
  #  slides = {./conf-presentations/presentation_IJCNN2009}
  # }
  # 
  # @article{chamroukhi_et_al_NN2009,
  # 	Address = {Oxford, UK, UK},
  # 	Author = {Chamroukhi, F. and Sam\'{e}, A. and Govaert, G. and Aknin, P.},
  # 	Date-Added = {2014-10-22 20:08:41 +0000},
  # 	Date-Modified = {2014-10-22 20:08:41 +0000},
  # 	Journal = {Neural Networks},
  # 	Number = {5-6},
  # 	Pages = {593--602},
  # 	Publisher = {Elsevier Science Ltd.},
  # 	Title = {Time series modeling by a regression approach based on a latent process},
  # 	Volume = {22},
  # 	Year = {2009},
  # 	url  = {https://chamroukhi.users.lmno.cnrs.fr/papers/Chamroukhi_Neural_Networks_2009.pdf}
  # 	}
  # 
  ##########################################################################################
  options(warn=-1)
  
  if (nargs()<9){verbose =0}
  if (nargs()<8){threshold = 1e-6}
  if (nargs()<7){max_iter_EM = 1500}
  if (nargs()<6){total_EM_tries = 1}
  if (nargs()<5){total_EM_tries = 1
  type_variance='hetereskedastic'}
  if (type_variance == 'homoskedastic'){homoskedastic =1}
  else if (type_variance == 'hetereskedastic'){homoskedastic=0}
  else{stop('The type of the model variance should be : "homoskedastic" or "hetereskedastic"')}
  
  #Chargement de toutes les fonctions ‡ utiliser
  source("designmatrix.R")
  source("init_hmmr.R")
  source("forwards_backwards.R")
  source("normalise.R")
  source("mk_stochastic.R")
  source("MAP.R")
  source("hmm_process.R")
  
  if (ncol(y)!=1){
    y=t(y)
  } 
  
  m = length(y)#length(y)

  X = designmatrix(x,p)#design matrix
  P = ncol(X)# here P is p+1
  I = diag(P)# define an identity matrix, in case of a Bayesian regularization for regression
  
  #
  best_loglik = -10000000
  nb_good_try=0
  total_nb_try=0
  cputime_total=c()
  
  while (nb_good_try < total_EM_tries){
    start_time = Sys.time()
    if (total_EM_tries>1){
      paste('EM try n∞',(nb_good_try+1))
    }
    total_nb_try=total_nb_try+1
    
    ## EM Initializaiton step
    ## Initialization of the Markov chain params, the regression coeffs, and the variance(s)
    HMMR = init_hmmr(X, y, K, type_variance, nb_good_try+1)
    
    # calculare the initial post probs  (tau_tk) and joint post probs (xi_ikl)
    
    #f_tk = hmmr.stats.f_tk; # observation component densities: f(yt|zt=k)
    prior = HMMR[[1]]
    trans_mat = HMMR[[3]]
    Mask = HMMR[[2]]
    betak = HMMR[[4]][[1]]
    if (homoskedastic==1){
      sigma2 = HMMR[[4]][[2]]
    }
    else {
      sigma2k = HMMR[[4]][[3]]
    }
    #
    iter = 0
    prev_loglik = -100000
    converged = 0
    top = 0
    
    #
    log_f_tk = matrix(c(0),m,K)
    muk = matrix(c(0),m,K)
    #
    ## EM
    stored_loglik=c()
    while ((iter <= max_iter_EM) & (converged!=1)){
      ## E step : calculate tge tau_tk (p(Zt=k|y1...ym;theta)) and xi t_kl (and the log-likelihood) by
      #  forwards backwards (computes the alpha_tk et beta_tk)
      
      # observation likelihoods
      for (k in 1:K){
        mk = X%*%betak[,k]
        muk[,k] = mk
        # the regressors means
        if (homoskedastic==1){
          sk = sigma2
        }
        else{
          sk = sigma2k[k]
        }
        z=((y - mk)^2)/sk
        log_f_tk[,k] =  -0.5*matrix(c(1),m,1)%*%(log(2*pi)+log(sk)) - 0.5*z#log(gaussienne)
        
      }
      for (k in 1:K){
        for (i in 1:nrow(log_f_tk)){
          log_f_tk[i,k]  = min(log_f_tk[i,k],log(.Machine$double.xmax))
          log_f_tk[i,k] = max(log_f_tk[i,k] ,log(.Machine$double.xmin))
        }
      }
      f_tk = exp(log_f_tk)
      
      #fprintf(1, 'forwards-backwards ');
      #[tau_tk, xi_tkl, alpha_tk, beta_tk, loglik] = 
      fb=forwards_backwards(prior, trans_mat , f_tk )
      tau_tk=fb[[1]]
      xi_tkl=fb$xi_tkl
      alpha_tk=fb[[3]]
      beta_tk=fb[[4]]
      loglik=fb[[5]]
      
      ## M step
      #  updates of the Markov chain parameters
      # initial states prob: P(Z_1 = k)
      prior = normalise(tau_tk[1,])[[1]]
      # transition matrix: P(Zt=i|Zt-1=j) (A_{k\ell})
      #print(cbind(apply(xi_tkl[,,1],2,sum),apply(xi_tkl[,,2],2,sum),apply(xi_tkl[,,3],2,sum)))
      # print(xi_tkl[,,1])

      trans_mat = round(mk_stochastic(cbind(apply(xi_tkl[,,1],2,sum),
                                      apply(xi_tkl[,,2],2,sum),
                                      apply(xi_tkl[,,3],2,sum),
                                      apply(xi_tkl[,,4],2,sum),
                                      apply(xi_tkl[,,5],2,sum))),4)
      
      # for segmental HMMR: p(z_t = k| z_{t-1} = \ell) = zero if k<\ell (no back) of if k >= \ell+2 (no jumps)
      trans_mat = mk_stochastic(Mask*trans_mat)
      ##  update of the regressors (reg coefficients betak and the variance(s) sigma2k)

      s = 0 # if homoskedastic
      for (k in 1:K){
        wieghts = tau_tk[,k]
        
        nk = sum(wieghts)# expected cardinal nbr of state k
        Xk = X*(sqrt(wieghts)%*%matrix(c(1),1,P))#[n*(p+1)]
        yk=y*(sqrt(wieghts))# dimension :[(nx1).*(nx1)] = [nx1]
        
        # reg coefficients
        lambda = 1e-9 # if a bayesian prior on the beta's
        bk = (solve(t(Xk)%*%Xk + lambda*diag(P))%*%t(Xk))%*%yk
        betak[,k] = bk
        
        # variance(s)
        z = sqrt(wieghts)*(y-X%*%bk)
        sk = t(z)%*%z
        if (homoskedastic==1){
          s = (s+sk)
          sigma2 = s/m
        }
        else{
          sigma2k[k] = sk/nk
        }
      }
      ## En of an EM iteration
      iter =  iter + 1
      
      # test of convergence
      loglik = loglik + log(lambda)
      
      if (verbose==1){
        paste('HMM_regression | EM   : Iteration :', iter,' Log-likelihood : ', loglik)
      }
      
      if (prev_loglik-loglik > 1e-4){
        top = top+1;
        if (top==10){
          stop(paste('!!!!! The loglikelihood is decreasing from',prev_loglik,' to ',loglik))
        }
      }
      converged = (abs(loglik - prev_loglik)/abs(prev_loglik) < threshold)
      stored_loglik[iter] <-loglik
      print(stored_loglik)
      prev_loglik = loglik
      end_time=Sys.time()
      cputime_total[nb_good_try+1]=c(end_time-start_time)
    }
    
    reg_param=list(betak=betak)
    
    if (homoskedastic==1){
      reg_param = append(reg_param,list(sigma2=sigma2))
    }
    else{
      reg_param = append(reg_param,list(sigma2k=sigma2k))
    }
    
    hmmr=list(prior=prior,trans_mat=trans_mat,reg_param=reg_param)
    
    # Estimated parameter vector (Pi,A,\theta)
    if (homoskedastic==1){
      parameter_vector=c(prior, trans_mat[Mask!=0],c(betak[1:P,]), sigma2)
      nu = K-1 + K*(K-1) + K*(p+1) + 1#length(parameter_vector);#
    }
    else{
      parameter_vector=c(prior, trans_mat[Mask!=0],c(betak[1:P,]), sigma2k)
      nu = K-1 + K*(K-1) + K*(p+1) + K #length(parameter_vector);#
    }
    
    stats=list(nu=nu,parameter_vector=parameter_vector,tau_tk=tau_tk,alpha_tk=alpha_tk,beta_tk=beta_tk,
               xi_tkl=xi_tkl,f_tk=f_tk,log_f_tk=log_f_tk,loglik=loglik,stored_loglik=stored_loglik,X=X)
    
    if (total_EM_tries>1){
      paste('loglik_max = ',loglik)
    }
    #
    #
    if (length(hmmr$reg_param$betak)!=0){
      nb_good_try=nb_good_try+1
      total_nb_try=0
      if (loglik > best_loglik){
        best_hmmr = hmmr
        best_loglik = loglik
      }
    }
    
    if (total_nb_try > 500){
      paste('can',"'",'t obtain the requested number of classes')
      hmmr=NULL
      return(hmmr)
    }
    
    
  }#End of the EM runs
  
  hmmr = best_hmmr
  
  #
  if (total_EM_tries>1){
    paste('best_loglik:  ',stats$loglik)
  }
  #
  #
  stats=append(stats,list(cputime=mean(cputime_total),cputime_total=cputime_total))
  
  ## Smoothing state sequences : argmax(smoothing probs), and corresponding binary allocations partition
  stats=append(stats,MAP(stats$tau_tk)) #[hmmr.stats.klas, hmmr.stats.Zik ]
  
  # #  compute the sequence with viterbi
  # #[path, ~] = viterbi_path(hmmr.prior, hmmr.trans_mat, hmmr.stats.fik');
  # #hmmr.stats.viterbi_path = path;
  # #hmmr.stats.klas = path;
  # ###################
  # 
  # # ## determination des temps de changements (les fonti√®tres entres les
  # # ## classes)
  # # nk=sum(hmmr.stats.Zik,1);
  # # for k = 1:K
  # #     tk(k) = sum(nk(1:k));
  # # end
  # # hmmr.stats.tk = [1 tk];
  # 
  # ## sate sequence prob p(z_1,...,z_n;\pi,A)
  state_probs = hmm_process(hmmr$prior, hmmr$trans_mat, m)
  stats = append(stats,list(state_probs=state_probs))
  
  ### BIC, AIC, ICL
  stats = append(stats,list(BIC=(stats$loglik - (stats$nu*log(m)/2))))
  stats = append(stats,list(AIC=(stats$loglik - stats$nu)))
  # # CL(theta) : Completed-data loglikelihood
  # sum_t_log_Pz_ftk = sum(hmmr.stats.Zik.*log(state_probs.*hmmr.stats.f_tk), 2);
  # comp_loglik = sum(sum_t_log_Pz_ftk(K:end));
  # hmmr.stats.comp_loglik = comp_loglik;
  # hmmr.stats.ICL = comp_loglik - (nu*log(m)/2);
  
  ## predicted, filtered, and smoothed time series
  stats = append(stats,list(regressors = round(X%*%hmmr$reg_param$betak,4)))
  
  # prediction probs   = Pr(z_t|y_1,...,y_{t-1})
  predict_prob = matrix(c(0),m,K)
  predict_prob[1,] = hmmr$prior#t=1 p (z_1)
  predict_prob[2:m,] = round((stats$alpha_tk[(1:(m-1)),]%*%hmmr$trans_mat)/(apply(stats$alpha_tk[(1:(m-1)),],1,sum)%*%matrix(c(1),1,K)),5)#t =2,...,n
  stats = append(stats,list(predict_prob = predict_prob))
  
  # predicted observations
  stats = append(stats,list(predicted = apply(round(stats$predict_prob*stats$regressors,5),1,sum)))#pond par les probas de prediction
  
  # filtering probs  = Pr(z_t|y_1,...,y_t)
  stats = append(stats,list(filter_prob = round(stats$alpha_tk/(apply(stats$alpha_tk,1,sum)%*%matrix(c(1),1,K)),5)))#normalize(alpha_tk,2);
  
  # filetered observations
  stats = append(stats,list(filtered = apply(round(stats$filter_prob*stats$regressors,5), 1,sum)))#pond par les probas de filtrage
  
  ### smoothed observations
  stats = append(stats,list(smoothed_regressors = stats$tau_tk*stats$regressors))
  stats =append(stats,list(smoothed = apply(stats$smoothed_regressors, 1,sum)))

  hmmr=c(hmmr,stats=list(stats))
  return(hmmr)
}

