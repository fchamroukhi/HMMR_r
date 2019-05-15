ModelHMMR <- setRefClass(
  "ModelHMMR",
  contains = "FData",
  # Define the fields
  fields = list(
    K = "numeric", # number of regimes
    p = "numeric", # dimension of beta (order of polynomial regression)
    variance_type = "numeric",
    nu = "numeric" # degree of freedom
  )
)

ModelHMMR <- function(fData, K, p, variance_type) {
  if (variance_type == variance_types$homoskedastic) {
    nu <<- K - 1 + K * (K - 1) + K * (p + 1) + 1
  }
  else{
    nu <<- K - 1 + K * (K - 1) + K * (p + 1) + K
  }

  new(
    "ModelHMMR",
    Y = fData$Y,
    X = fData$X,
    m = fData$m,
    n = fData$n,
    K = K,
    p = p,
    variance_type = variance_type,
    nu = nu
  )
}
