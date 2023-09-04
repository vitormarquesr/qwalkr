#' The All-Ones Matrix
#'
#' @param n the order of the matrix.
#'
#' @returns A square matrix of order \eqn{n} in which every entry is
#'   equal to 1. The all-ones matrix is given by
#'   \eqn{J_{n\ x\ n} = 1_{n\ x\ 1}1_{n\ x\ 1}^T}.
#'
#' @export
#'
#' @examples
#' # Return the all-ones matrix of order 5.
#' J(5)
J <- function(n){
  return (matrix(1, nrow=n, ncol=n))
}

#' The Trace of a Matrix
#'
#' @param A a square matrix.
#'
#' @returns If \eqn{A} has order \eqn{n},
#'   then \eqn{tr(A) = \sum_{i=1}^{n}a_{ii}}.
#'
#' @export
#'
#' @examples
#' A <- rbind(1:5, 2:6, 3:7)
#'
#' # Calculate the trace of A
#' tr(A)
tr <- function(A){
  return (sum(diag(A)))
}

#' The Trace Inner Product of Matrices
#'
#' @param A,B square matrices.
#'
#' @returns The trace inner product on \eqn{Mat_{n\ x\ n}(\mathbb{C})} is
#'   defined as
#'
#'   \deqn{\langle A, B \rangle := tr(A^*B)}
#' @export
#'
#' @examples
#' A <- rbind(1:5, 2:6, 3:7)
#' B <- rbind(7:11, 8:12, 9:13)
#'
#' # Calculate the trace inner product of A and B
#' trdot(A, B)
trdot <- function(A, B){
  return (tr(Conj(t(A)) %*% B))
}



