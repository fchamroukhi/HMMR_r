source("R/enums.R")
source("R/mk_stochastic.R")
# source("R/utils.R")
# source("R/IRLS.R")

ParamHMMR <- setRefClass(
  "ParamHMMR",
  fields = list(prior = "matrix",
                trans_mat = "matrix",
                beta = "matrix",
                sigma = "matrix",
                mask = "matrix"),
  methods = list(

    init_hmmr = function(modelHMMR, phi, try_algo = 1) {
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

      ## Tnitialisation en tenant compte de la contrainte:

      # Initialisation de la matrice des transitions
      Mask <- 0.5 * diag(modelHMMR$K) # masque d'ordre 1

      for (k in 1:(modelHMMR$K-1)) {
        ind <- which(Mask[k,] != 0)
        Mask[k,ind+1] <- 0.5
      }
      trans_mat <<- Mask
      mask <<- Mask

      # Initialisation de la loi initiale de la variable cachee
      prior <<- matrix(c(1,rep(0,modelHMMR$K-1)))

      #  Initialisation des coeffecients de regression et des variances.
      init_hmmr_regressors(phi, modelHMMR, try_algo)

    },

    init_hmmr_regressors = function (phi, modelHMMR, try_algo = 1){

      if(try_algo ==1){# uniform segmentation into K contiguous segments, and then a regression
        zi = round(modelHMMR$m/modelHMMR$K)-1

        s=0
        for (k in 1:modelHMMR$K){
          yk <- modelHMMR$Y[((k-1)*zi+1):(k*zi)]
          Xk <- phi[((k-1)*zi+1):(k*zi),]

          beta[,k] <<- solve(t(Xk)%*%Xk +(10^-4)*diag(modelHMMR$p + 1))%*%t(Xk)%*%yk #regress(yk,Xk); # for a use in octave, where regress doesnt exist
          muk <- Xk%*%beta[,k]
          sk <- t(yk-muk)%*%(yk-muk)
          if (modelHMMR$variance_type==variance_types$homoskedastic){
            s<- (s+sk)
            sigma <<- s/modelHMMR$m
          }
          else {
            sigma[k] <<- sk/length(yk)
          }
        }
      }
      else{    # random segmentation into contiguous segments, and then a regression
        Lmin <- modelHMMR$p + 1 + 1#minimum length of a segment #10
        tk_init <- rep(0,modelHMMR$K)
        tk_init <- t(tk_init)
        tk_init[1] <- 0
        K_1<-modelHMMR$K
        for (k in 2:modelHMMR$K){
          K_1 <- K_1-1
          temp <- seq(tk_init[k-1]+Lmin,modelHMMR$m-K_1*Lmin)
          ind <- sample(1:length(temp),length(temp))
          tk_init[k]<- temp[ind[1]]
        }
        tk_init[K+1]<-modelHMMR$m

        s=0#
        for (k in 1:modelHMMR$K){
          i <- tk_init[k]+1
          j <- tk_init[k+1]
          yk <- modelHMMR$Y[i:j]
          Xk <- phi[i:j,]
          beta[,k] <<- solve(t(Xk)%*%Xk + 1e-4*diag(modelHMMR$p + 1))%*%t(Xk)%*%yk #regress(yk,Xk); # for a use in octave, where regress doesnt exist
          muk <- Xk%*%beta[,k]
          sk <- t(yk-muk)%*%(yk-muk)

          if (modelHMMR$variance_type==variance_types$homoskedastic){
            s <- s+sk
            sigma <<- s/modelHMMR$m

          }
          else{
            sigma[k] <<- sk/length(yk)
          }
        }
      }
    },

    MStep = function(modelHMMR, statHMMR, phi) {

      #  updates of the Markov chain parameters
      # initial states prob: P(Z_1 = k)
      prior <<- matrix(normalize(statHMMR$tau_tk[1,])$M)
      # transition matrix: P(Zt=i|Zt-1=j) (A_{k\ell})
      #print(cbind(apply(xi_tkl[,,1],2,sum),apply(xi_tkl[,,2],2,sum),apply(xi_tkl[,,3],2,sum)))
      # print(xi_tkl[,,1])
      trans = NULL
      for (k in 1:modelHMMR$K){
        trans = cbind(trans,apply(statHMMR$xi_tkl[,,k],2,sum))
      }

      trans_mat <<- round(mk_stochastic(trans), 4)

      # for segmental HMMR: p(z_t = k| z_{t-1} = \ell) = zero if k<\ell (no back) of if k >= \ell+2 (no jumps)
      trans_mat <<- mk_stochastic(mask * trans_mat)
      ##  update of the regressors (reg coefficients betak and the variance(s) sigma2k)

      s = 0 # if homoskedastic
      for (k in 1:modelHMMR$K){
        weights = statHMMR$tau_tk[,k]

        nk = sum(weights)# expected cardinal nbr of state k
        Xk = phi*(sqrt(weights)%*%matrix(1,1,modelHMMR$p + 1))#[n*(p+1)]
        yk=modelHMMR$Y*(sqrt(weights))# dimension :[(nx1).*(nx1)] = [nx1]

        # reg coefficients
        lambda = 1e-9 # if a bayesian prior on the beta's
        bk = (solve(t(Xk)%*%Xk + lambda*diag(modelHMMR$p + 1))%*%t(Xk))%*%yk
        beta[,k] <<- bk

        # variance(s)
        z = sqrt(weights)*(modelHMMR$Y-phi%*%bk)
        sk = t(z)%*%z
        if (modelHMMR$variance_type == variance_types$homoskedastic){
          s = (s+sk)
          sigma <<- s/modelHMMR$m
        }
        else{
          sigma[k] <<- sk/nk
        }
      }

    }
  )
)

ParamHMMR <- function(modelHMMR) {
  prior <- matrix(NA, ncol = modelHMMR$K - 1)
  trans_mat <- matrix(NA, modelHMMR$K, modelHMMR$K)
  beta <- matrix(NA, modelHMMR$p + 1, modelHMMR$K)
  if (modelHMMR$variance_type == variance_types$homoskedastic) {
    sigma <- matrix(NA)
  }
  else{
    sigma <- matrix(NA, modelHMMR$K)
  }
  mask <- matrix(NA, modelHMMR$K, modelHMMR$K)
  new("ParamHMMR", prior = prior, trans_mat = trans_mat, beta = beta, sigma = sigma, mask = mask)
}
