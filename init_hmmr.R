##################################################################################################################################################
init_hmmr <- function(X, y, K, type_variance, EM_try){
# function hmmr =  init_hmmr(X, y, K, type_variance, EM_try)
# init_hmmr initialize the parameters of a Gaussian Hidden Markov Model
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
#         and some stats: like the Mask for a segmental model
#
#
# Faicel Chamroukhi, first version in November 2008
##################################################################################
  if (type_variance == 'homoskedastic'){
    homoskedastic =1
    }
  else{
    homoskedastic=0
    }
  m = length(y)
## Tnitialisation en tenant compte de la contrainte:

# Initialisation de la matrice des transitions
  Mask = (0.5)*diag(K) #masque d'ordre 1
  for (k in 1:(K-1)) {
    ind = which(Mask[k,] != 0)
    Mask[k,ind+1] = 0.5
  }
  hmmr.trans_mat = Mask

# Initialisation de la loi initiale de la variable cachee
  hmmr.prior = t(c(1,rep(0,K-1)))
  hmmr.stats.Mask = Mask

#  Initialisation des coeffecients de regression et des variances.
  reg_param = init_hmmr_regressors(X, y, K, type_variance, EM_try)
  
  return(list(hmmr.prior,hmmr.stats.Mask,hmmr.trans_mat,reg_param))
}
################################################################################
init_hmmr_regressors <- function (X, y, K, type_variance, EM_try){
  if (type_variance == 'homoskedastic'){
    homoskedastic =1
    }
  else{
     homoskedastic = 0
  }
  m= nrow(X)
  P= ncol(X)
   #m = length(y);
  betak=matrix(nrow=P,ncol=K)
  sigma2k=matrix(nrow=1,ncol=K)
  sigma2=c(0)
   
  if(EM_try ==1){# uniform segmentation into K contiguous segments, and then a regression
     zi = round(m/K)-1
     
     s=0
     for (k in 1:K){
       yk = y[((k-1)*zi+1):(k*zi)]
       Xk = X[((k-1)*zi+1):(k*zi),]
       
       betak[,k] = solve(t(Xk)%*%Xk +(10^-4)*diag(P))%*%t(Xk)%*%yk #regress(yk,Xk); # for a use in octave, where regress doesnt exist
       muk = Xk%*%betak[,k]
       sk = t(yk-muk)%*%(yk-muk)
       if (homoskedastic==1){
         s= (s+sk)
         sigma2 = s/m
       }
       else {
         sigma2k[k] = sk/length(yk)
       }
     }
   }
  else{    # random segmentation into contiguous segments, and then a regression
    Lmin= P+1#minimum length of a segment #10
    tk_init = rep(c(0),K)
    tk_init = t(tk_init)
    tk_init(1) = 0
    K_1=K
    for (k in 2:K){
      K_1 = K_1-1
      temp = tk_init(k-1)+sum(Lmin:m)-K_1*Lmin
      ind = randperm(length(temp))
      tk_init[k]= temp(ind(1))
    }
    tk_init[k+1] = m
    
    s=0#
    for (k in 1:K){
      i = tk_init[k]+1
      j = tk_init[k+1]
      yk = y(i:j)#y((k-1)*zi+1:k*zi);
      Xk = X[i:j,]#X((k-1)*zi+1:k*zi,:);
      betak[,k] = solve(t(Xk)%*%Xk + 1e-4*diag(P))%*%t(Xk)%*%yk #regress(yk,Xk); # for a use in octave, where regress doesnt exist
      muk = Xk%*%betak[,k]
      sk = t(yk-muk)*(yk-muk)
      
      if (homoskedastic==1){
        s = s+sk
        sigma2 = s/m
      }
      else{
        sigma2k[k] = sk/length(yk)
      }
    }
  }
  return(list(betak=betak,sigma=sigma2,sigma2k=sigma2k))
}

