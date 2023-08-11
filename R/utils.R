
equal <- function(x,y, tol = .Machine$double.eps^0.5){
  return (abs(x-y) < tol)
}

unique_eigvals <- function(x, ...){
  n <- length(x)
  return (!equal(x, c(Inf, x[seq_len(n-1)]), ...))
}

mult_eigvals <- function(tags){
  n <- length(tags)
  return (diff(c(which(tags), n+1)))
}

index_eigspace <- function(mult, id){
  idx_start <- sum(mult[seq_len(id-1)]) + 1
  return (idx_start:(idx_start + mult[id]-1))
}

out_of_bounds <- function(id, lower, upper){
  return (any((id < lower) | (id > upper)))
}


