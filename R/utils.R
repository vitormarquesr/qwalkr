
equal <- function(x,y, tol = .Machine$double.eps^0.5){
  return(abs(x-y) < tol)
}

label_unique <- function(x, tol = .Machine$double.eps^0.5){
  n <- length(x)
  labeled <- vector(mode="integer", length=n)
  idx <- 0L

  for (i in 1:n){
    if ((i == 1) || !equal(x[i], x[i-1L], tol=tol)){
      idx <- idx + 1L
    }
    labeled[[i]] <- idx
  }

  return(labeled)
}

index_eigspace <- function(mult, id){
  idx_start <- sum(mult[seq_len(id-1)]) + 1
  return(idx_start:(idx_start + mult[id]-1))
}

out_of_bounds <- function(id, lower, upper){
  return (any((id < lower) | (id > upper)))
}


