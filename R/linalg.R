#' Spectral Decomposition of a Hermitian Matrix
#'
#' @description `spectral()` is a wrapper around [base::eigen()] designed for Hermitian matrices,
#'   which can handle repeated eigenvalues.
#'
#' @param S a Hermitian matrix. *Obs*: The matrix is always assumed to be Hermitian,
#'   and only its lower triangle (diagonal included) is used.
#' @param multiplicity if `TRUE` (default), tries to infer eigenvalue multiplicity. If set to
#'   `FALSE`, each eigenvalue is considered unique with multiplicity one.
#' @param tol two eigenvalues `x`, `y` are considered equal if `abs(x-y) < tol`. Defaults to
#'   `tol=.Machine$double.eps^0.5`.
#' @param ... further arguments passed on to [base::eigen()]
#'
#' @returns The spectral decomposition of `S` is returned as a list with components
#'   \item{eigvals}{vector containing the unique eigenvalues of `S` in *decreasing* order.}
#'   \item{multiplicity}{multiplicities of the eigenvalues in `eigvals`.}
#'   \item{eigvectors}{a `nrow(S) x nrow(S)` unitary matrix whose columns are eigenvectors ordered
#'   according to `eigvals`. Note that there may be more eigenvectors than eigenvalues if
#'   `multiplicity=TRUE`, however eigenvectors of the same eigenspace are next to each other.}
#'
#'   The Spectral Theorem ensures the eigenvalues of `S` are real and that the vector space
#'   admits an orthonormal basis consisting of eigenvectors of `S`. Thus, if `s <- spectral(S)`,
#'   and `V <- s$eigvectors; lam <- s$eigvals`, then
#'
#'   \deqn{S = V \Lambda V^{*}}
#'
#'   where \eqn{\Lambda =\ }`diag(rep(lam, times=s$multiplicity))`
#'
#' @seealso [base::eigen()], [qwalkr::extractEIGSPACE.spectral],
#'   [qwalkr::extractPROJ.spectral], [qwalkr::extractSCHUR.spectral]
#'
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
spectral <- function(S, multiplicity = TRUE, tol=.Machine$double.eps^0.5, ...){

  decomp <- eigen(S, symmetric=TRUE, ...)

  n <- nrow(S)

  idx_eigvals <- if (multiplicity) unique_eigvals(decomp$values, tol=tol) else rep(TRUE, n)

  output <- list(eigvals = decomp$values[idx_eigvals],
                 multiplicity = mult_eigvals(idx_eigvals),
                 eigvectors = decomp$vectors)

  class(output) <- "spectral"

  return (output)
}

#' Generic S3 method extractEIGSPACE
#'
#' @param object an object containing the spectral decomposition of a matrix.
#' @param ... further arguments passed to or from other methods.
#'
#' @returns The eigenbasis of the desired eigenspace.
#' @seealso [qwalkr::extractPROJ], [qwalkr::extractSCHUR],
#'   [qwalkr::extractEIGSPACE.spectral]
#' @export
#'
extractEIGSPACE <- function(object, ...) UseMethod("extractEIGSPACE")


#' extractEIGSPACE method for spectral objects
#'
#' @param object an object of class spectral.
#' @param id index for the desired eigenspace according to the ordered (decreasing) spectra.
#' @param ... further arguments passed to or from other methods.
#' @returns A matrix whose columns form the orthonormal eigenbasis.
#'
#'   If `s <- spectral(A)` and `V <- s$eigvectors`, then the extracted eigenspace
#'   \eqn{V_{id}} is some submatrix `V[, _]`.
#'
#' @export
#' @seealso [qwalkr::spectral], [qwalkr::extractEIGSPACE]
#'
#' @examples
#' # Spectra is {2, -1} with multiplicities one and two respectively.
#' decomp <- spectral(matrix(c(0,1,1,1,0,1,1,1,0), nrow=3))
#'
#' # Returns the two orthonormal eigenvectors corresponding to the eigenvalue -1.
#' extractEIGSPACE(decomp, id=2)
#'
#' # Returns the eigenvector corresponding to the eigenvalue 2.
#' extractEIGSPACE(decomp, id=1)
#'
extractEIGSPACE.spectral <- function(object, id, ...){
  if (out_of_bounds(id, 1, length(object$multiplicity))){
    warning("Index out of bounds! Check the length of the spectra.")
    return(NULL)
  }

  return (object$eigvectors[, index_eigspace(object$multiplicity, id), drop=FALSE])
}

#' Generic S3 method extractPROJ
#'
#' @param object an object containing the spectral decomposition of a matrix.
#' @param ... further arguments passed to or from other methods.
#'
#' @returns The orthogonal projector of the desired eigenspace.
#' @seealso [qwalkr::extractEIGSPACE], [qwalkr::extractSCHUR],
#'    [qwalkr::extractPROJ.spectral]
#' @export
#'
extractPROJ <- function(object, ...) UseMethod("extractPROJ")

#' extractPROJ method for spectral objects
#'
#' @param object an object of class spectral.
#' @param id index for the desired eigenspace according to the ordered (decreasing) spectra.
#' @param ... further arguments passed to or from other methods.
#'
#' @returns The orthogonal projector of the desired eigenspace.
#'
#'   A Hermitian matrix `S` admits the spectral decomposition \eqn{S = \sum_{r}\lambda_r E_r}
#'   such that \eqn{E_r} is the orthogonal projector onto the \eqn{\lambda_r}-eigenspace. If \eqn{V_{id}}
#'   is the matrix associated to the eigenspace, then
#'
#'   \deqn{E_{id} = V_{id}V_{id}^*}
#'
#' @seealso [qwalkr::spectral], [qwalkr::extractPROJ]
#'
#' @export
#'
#' @examples
#' # Spectra is {2, -1} with multiplicities one and two respectively.
#' decomp <- spectral(matrix(c(0,1,1,1,0,1,1,1,0), nrow=3))
#'
#' # Returns the projector associated to the eigenvalue -1.
#' extractPROJ(decomp, id=2)
#'
#' # Returns the projector associated to the eigenvalue 2.
#' extractPROJ(decomp, id=1)
#'
extractPROJ.spectral <- function(object, id, ...){
  if (out_of_bounds(id, 1, length(object$multiplicity))){
    warning("Index out of bounds! Check the length of the spectra.")
    return(NULL)
  }

  A <- extractEIGSPACE.spectral(object, id)
  return (A %*% Conj(t(A)))
}

#' Generic S3 method extractSCHUR
#'
#' @param object an object containing the spectral decomposition of a matrix.
#' @param ... further arguments passed to or from other methods.
#'
#' @returns The Schur product of eigenprojectors.
#' @seealso [qwalkr::extractEIGSPACE], [qwalkr::extractSCHUR],
#'    [qwalkr::extractSCHUR.spectral]
#' @export
#'
extractSCHUR <- function(object, ...) UseMethod("extractSCHUR")

#' extractSCHUR method for spectral objects
#'
#' @param object an object of class spectral.
#' @param id1 index for the first eigenspace according to the ordered (decreasing) spectra.
#' @param id2 index for the second eigenspace according to the ordered (decreasing) spectra. If not provided,
#'    it takes the same value as `id1`.
#' @param ... further arguments passed to or from other methods.
#'
#' @returns The Schur product of the corresponding eigenprojectors, \eqn{E_{id_1} \circ E_{id_2}}.
#'
#' @seealso [qwalkr::spectral], [qwalkr::extractSCHUR]
#'
#' @export
#'
#' @examples
#' # Spectra is {2, -1} with multiplicities one and two respectively.
#' decomp <- spectral(matrix(c(0,1,1,1,0,1,1,1,0), nrow=3))
#'
#' # Returns the Schur product between the 2-projector and -1-projector.
#' extractSCHUR(decomp, id1=2, id2=1)
#'
#' # Returns the Schur square of the 2-projector.
#' extractSCHUR(decomp, id1=1, id2=1)
#'
#' # Also returns the Schur square of the 2-projector
#' extractSCHUR(decomp, id1=1)
#'
extractSCHUR.spectral <- function(object, id1, id2=NULL, ...){

  id2 <- if(is.null(id2)) id1 else id2

  if (out_of_bounds(c(id1, id2), 1, length(object$multiplicity))){
    warning("Index out of bounds! Check the length of the spectra.")
    return(NULL)
  }
  E_r <- extractPROJ.spectral(object, id1)
  E_s <- extractPROJ.spectral(object, id2)
  return (E_r * E_s)
}



