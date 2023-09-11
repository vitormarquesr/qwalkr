#' The All-Ones Matrix
#'
#' @description
#' Returns the all-ones matrix of order `n`.
#'
#' @param n the order of the matrix.
#'
#' @returns A square matrix of order \eqn{n} in which every entry is
#'   equal to 1. The all-ones matrix is given by
#'   \eqn{J_{n\ x\ n} = 1_{n\ x\ 1}1_{n\ x\ 1}^T}.
#'
#' @export
#'
#' @seealso [qwalkr::tr()], [qwalkr::trdot()], [qwalkr::cartesian()]
#'
#' @examples
#' # Return the all-ones matrix of order 5.
#' J(5)
J <- function(n){
  return (matrix(1, nrow=n, ncol=n))
}

#' The Trace of a Matrix
#'
#' @description
#' Computes the trace of a matrix `A`.
#'
#' @param A a square matrix.
#'
#' @returns If \eqn{A} has order \eqn{n},
#'   then \eqn{tr(A) = \sum_{i=1}^{n}a_{ii}}.
#'
#' @export
#' @seealso [qwalkr::J()], [qwalkr::trdot()], [qwalkr::cartesian()]
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
#' @description
#' Computes the trace inner product of two matrices `A` and `B`.
#'
#' @param A,B square matrices.
#'
#' @returns The trace inner product on \eqn{Mat_{n\ x\ n}(\mathbb{C})} is
#'   defined as
#'
#'   \deqn{\langle A, B \rangle := tr(A^*B)}
#' @export
#' @seealso [qwalkr::J()], [qwalkr::tr()], [qwalkr::cartesian()]
#'
#' @examples
#' A <- rbind(1:5, 2:6, 3:7)
#' B <- rbind(7:11, 8:12, 9:13)
#'
#' # Compute the trace inner product of A and B
#' trdot(A, B)
trdot <- function(A, B){
  return (tr(Conj(t(A)) %*% B))
}

#' Adjacency Matrix of the Cartesian Product
#'
#' @description
#' Returns the adjacency matrix of the cartesian product of two graphs
#' given the adjacency matrix of each one, \eqn{G} and \eqn{H}.
#'
#' @param G adjacency matrix of the first graph.
#' @param H adjacency matrix of the second graph. If not provided,
#'   it takes the same value as `G`.
#'
#' @returns Let \eqn{A(G),\ A(H)} be the adjacency matrices of
#'   the graphs \eqn{G,\ H} such that \eqn{|V(G)| = n} and \eqn{|V(H)| = m},
#'   then the adjacency matrix of the cartesian product \eqn{G \times H} is
#'   given by
#'
#'   \deqn{A(G \times H) = A(G) \otimes I_{m\ x\ m} + I_{n\ x\ n} \otimes A(H)}
#'
#' @export
#' @seealso [qwalkr::J()], [qwalkr::tr()], [qwalkr::trdot()]
#'
#' @examples
#' P3 <- matrix(c(0,1,0,1,0,1,0,1,0), nrow=3)
#' K3 <- matrix(c(0,1,1,1,0,1,1,1,0), nrow=3)
#'
#' # Return the adjacency matrix of P3 X K3
#' cartesian(P3, K3)
#'
#' # Return the adjacency matrix of P3 X P3
#' cartesian(P3)
cartesian <- function(G, H=NULL){
  n1 <- nrow(G)
  if (is.null(H)){
    H <- G
  }
  n2 <- nrow(H)

  W <- kronecker(G, diag(n2)) + kronecker(diag(n1), H)

  return (W)
}


