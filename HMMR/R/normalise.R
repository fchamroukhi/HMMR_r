normalise <- function(M){
# NORMALISE Make the entries of a (multidimensional) array sum to 1
# [M, total] = normalise(M)
#
# This function uses the British spelling to avoid any confusion with
# 'normalize_pot', which is a method for potentials.
########################################################################
  total = sum(M)
  # Set any zeros to one before dividing
  d = total + as.numeric(total==0)
  M = M / d
  return(list(M = M, total = total))
}
