#' Spectral Decomposition of a Matrix
#'
#' @description `spectral()` is a wrapper around [base::eigen()] that can handle repeated eigenvalues.
#'
#' @param A a square matrix.
#' @param multiplicity if `TRUE` (default), tries to infer eigenvalue multiplicity. If set to
#'    `FALSE`, each eigenvalue is considered unique with multiplicity one.
#' @param tol two eigenvalues `x`, `y` are considered equal if `abs(x-y) < tol`. Defaults to
#'    `tol=.Machine$double.eps^0.5`.
#' @param ... further arguments passed on to [base::eigen()]
#'
#' @returns The spectral decomposition of `A` is returned as a list with components
#' \item{eigvals}{vector containing the unique eigenvalues of `A` in *decreasing* order.}
#' \item{multiplicity}{multiplicities of the eigenvalues in `eigvals`.}
#' \item{eigvectors}{a `nrow(A) x nrow(A)` matrix whose columns are eigenvectors ordered
#' according to `eigvals`. Note that there may be more eigenvectors than eigenvalues if
#' `multiplicity=TRUE`, however eigenvectors of the same eigenspace are next to each other.}
#'
#' @seealso [base::eigen()]
#' @export
#'
#' @examples
#' spectral(matrix(c(0,1,0,1,0,1,0,1,0), nrow=3))
#'
#' # Use "tol" to set the tolerance for numerical equality
#' spectral(matrix(c(0,1,0,1,0,1,0,1,0), nrow=3), tol=10e-5)
#'
#' # Use "multiplicity=FALSE" to force each eigenvalue to be considered unique
#' spectral(matrix(c(0,1,0,1,0,1,0,1,0), nrow=3), multiplicity = FALSE)
#'
#' # Use ... to pass additional arguments to eigen()
#' spectral(matrix(c(0,1,0,1,0,1,0,1,0), nrow=3), symmetric=TRUE)
#'
spectral <- function(A, multiplicity = TRUE, tol=.Machine$double.eps^0.5, ...){

  eigen_decomp <- eigen(A, ...)

  n <- nrow(A)


  output <- list(eigvals = eigen_decomp$values,
                 multiplicity = rep(1, n),
                 eigvectors = eigen_decomp$vectors)

  if (multiplicity){
    label_eigvals <- label_unique(output$eigvals, tol)

    output$eigvals <- output$eigvals[!duplicated(label_eigvals)]

    output$multiplicity <-  tabulate(label_eigvals)
  }

  class(output) <- "spectral"

  return (output)
}



