# User-freindly and flexible model anf algorithm for time series segmentation with a Regression
# model with a Hidden Markov Model Regression (HMMR).
# 
# 
# Hidden Markov Model Regression (HMMR) for segmentation of time series
# with regime changes. The model assumes that the time series is
# governed by a sequence of hidden discrete regimes/states, where each
# regime/state has Gaussian regressors as observations.
# The model parameters are estimated by MLE via the EM algorithm
# 
# Faicel Chamroukhi
# 
# Please cite the following papers for this code:
# 
# 
# @article{Chamroukhi-FDA-2018,
# 	Journal = {Wiley Interdisciplinary Reviews: Data Mining and Knowledge Discovery},
#	  Author = {Faicel Chamroukhi and Hien D. Nguyen},
#	  Note = {DOI: 10.1002/widm.1298.},
#	  Volume = {},
# 	Title = {Model-Based Clustering and Classification of Functional Data},
#	  Year = {2018},
#	  Month = {Dec},
#	  url =  {https://chamroukhi.com/papers/MBCC-FDA.pdf}
#	}
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
#   url = {https://chamroukhi.com/papers/chamroukhi_ijcnn2009.pdf},
#   slides = {./conf-presentations/presentation_IJCNN2009}
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
# 
# 
# 
# Faicel Chamroukhi Septembre 2008.
##############################################################################################
shell("cls") #efface la console
rm(list=ls())#efface les variables
#setwd("") #Repertoire de travail

#model specification
K = 5 # nomber of regimes (states)
p = 2 # dimension of beta' (order of the polynomial regressors)

# options
#type_variance = 'homoskedastic'
type_variance = 'hetereskedastic'
nbr_EM_tries = 1
max_iter_EM = 1500
threshold = 1e-6
verbose = 1
type_algo = 'EM'
#type_algo = 'CEM'
#type_algo = 'SEM'

## toy time series with regime changes
# y =[randn(100,1); 7+randn(120,1);4+randn(200,1); -2+randn(100,1); 3.5+randn(150,1);]'
# n = length(y)
# x = linspace(0,1,n)

library(R.matlab)
tmp = readMat("simulated_time_series.mat")
#tmp = readMat("real_time_series_1.mat")
#tmp = readMat("real_time_series_2.mat")
x = tmp$x
y = tmp$y

source("learn_hmmr.R")

# #model selection
# classe=2:7
# poly=2:3
# bic=matrix(c(0),length(classe),length(poly))
# current_BIC = -99999
# for (K in classe){
#   for (p in poly){
#     print(paste("K=",K,"p=",p))
#     hmmr = learn_hmmr(x, y, K, p, type_variance, nbr_EM_tries, max_iter_EM, threshold, verbose)
# 
#     if (hmmr$stats$BIC > current_BIC){
#       best_hmmr = hmmr
#       current_BIC = hmmr$stats$BIC
#     }
#     bic[(K-1),(p-1)] = hmmr$stats$BIC
#   }
# }
# 
# colors = c("red","blue","green","pink","cadetblue2","orange","blue4","chartreuse4","brown2","cadetblue4")
# x11()
# plot(classe,bic[,1],type="l",col=colors[1],ylab="BIC",xlab="K",main="S?lection de mod?le")
# for (p in 2:length(poly)){
#   lines(classe,bic[,p],type="l",col=colors[p])
# }
# legend(2,y=-1150,legend=c(min(poly):max(poly)),col=colors[1:length(poly)], lty=1:2, cex=0.8)

hmmr = learn_hmmr(x, y, K, p, type_variance, nbr_EM_tries, max_iter_EM, threshold, verbose)

source("show_HMMR_results.R")
# #yaxislim = c(240,600)
show_HMMR_results(x,y,hmmr)

# sample an HMMR
source("sample_hmmr.R")
#sample = sample_hmmr(x, hmmr$prior, hmmr$trans_mat, hmmr$reg_param$betak,hmmr$reg_param$sigma2k)

