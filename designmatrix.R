designmatrix <- function(x,p){
# constructs the desing matrix of a polynomial regression with degree p
#
#
# Faicel Chamroukhi
#######################################################################
  if (NCOL(x) != 1){ x=t(x)} #a column vector
  
  n = length(x)
  X=matrix(data=0, nrow=n, ncol=p+1)
  for (i in 0:p){
  X[,(i+1)] = x^(i)} # X = [1 x x^2 x^3 x^p]
  
  return(X)
}

