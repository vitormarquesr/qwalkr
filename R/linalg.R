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
#' @seealso [base::eigen()], [qwalkr::get_eigspace.spectral()],
#'   [qwalkr::get_eigproj.spectral()], [qwalkr::get_eigschur.spectral()],
#'   [qwalkr::act_eigfun.spectral()]
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

#' Extract an Eigenspace from an Operator
#'
#' @param object a representation of the operator.
#' @param ... further arguments passed to or from other methods.
#'
#' @returns A representation of the requested eigenspace.
#' @seealso [qwalkr::get_eigproj()], [qwalkr::get_eigschur()],
#'   [qwalkr::get_eigspace.spectral()]
#' @export
#' @examples
#' s <- spectral(rbind(c(0.5, 0.3), c(0.3,0.7)))
#'
#' get_eigspace(s, 1) #-> get_eigspace.spectral(...)
#'
get_eigspace <- function(object, ...) UseMethod("get_eigspace")


#' Extract an Eigenspace from a Hermitian Matrix
#'
#' @description
#' Get the eigenbasis associated with an eigenvalue based on the representation
#' of a Hermitian Matrix given by class `spectral`.
#'
#' @param object an instance of class `spectral`.
#' @param id index for the desired eigenspace according to the ordered (decreasing) spectra.
#' @param ... further arguments passed to or from other methods.
#' @returns A matrix whose columns form the orthonormal eigenbasis.
#'
#'   If `s <- spectral(A)` and `V <- s$eigvectors`, then the extracted eigenspace
#'   \eqn{V_{id}} is some submatrix `V[, _]`.
#'
#' @export
#' @seealso [qwalkr::spectral()], [qwalkr::get_eigspace()]
#'
#' @examples
#' # Spectra is {2, -1} with multiplicities one and two respectively.
#' decomp <- spectral(matrix(c(0,1,1,1,0,1,1,1,0), nrow=3))
#'
#' # Returns the two orthonormal eigenvectors corresponding to the eigenvalue -1.
#' get_eigspace(decomp, id=2)
#'
#' # Returns the eigenvector corresponding to the eigenvalue 2.
#' get_eigspace(decomp, id=1)
#'
get_eigspace.spectral <- function(object, id, ...){
  if (out_of_bounds(id, 1, length(object$multiplicity))){
    stop("Index out of bounds! Check the length of the spectra.")
  }

  return (object$eigvectors[, index_eigspace(object$multiplicity, id), drop=FALSE])
}

#' Extract an Eigen-Projector from an operator
#'
#' @param object a representation of the operator.
#' @param ... further arguments passed to or from other methods.
#'
#' @returns A representation of the requested eigen-projector.
#' @seealso [qwalkr::get_eigspace()], [qwalkr::get_eigschur()],
#'    [qwalkr::get_eigproj.spectral()]
#' @export
#' @examples
#' s <- spectral(rbind(c(0.5, 0.3), c(0.3,0.7)))
#'
#' get_eigproj(s, 1) #-> get_eigproj.spectral(...)
#'
get_eigproj <- function(object, ...) UseMethod("get_eigproj")

#' Extract an Eigen-Projector from a Hermitian Matrix
#'
#' @description
#' Get the orthogonal projector associated with an eigenspace based on the representation
#' of a Hermitian Matrix given by class `spectral`.
#'
#' @param object an instance of class `spectral`.
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
#' @seealso [qwalkr::spectral()], [qwalkr::get_eigproj()]
#'
#' @export
#'
#' @examples
#' # Spectra is {2, -1} with multiplicities one and two respectively.
#' decomp <- spectral(matrix(c(0,1,1,1,0,1,1,1,0), nrow=3))
#'
#' # Returns the projector associated to the eigenvalue -1.
#' get_eigproj(decomp, id=2)
#'
#' # Returns the projector associated to the eigenvalue 2.
#' get_eigproj(decomp, id=1)
#'
get_eigproj.spectral <- function(object, id, ...){
  if (out_of_bounds(id, 1, length(object$multiplicity))){
    stop("Index out of bounds! Check the length of the spectra.")
  }

  A <- get_eigspace.spectral(object, id)
  return (A %*% Conj(t(A)))
}

#' Extract a Schur Cross-Product from an Operator
#'
#' @param object a representation of the operator.
#' @param ... further arguments passed to or from other methods.
#'
#' @returns A representation of the requested Schur cross-product.
#' @seealso [qwalkr::get_eigspace()], [qwalkr::get_eigproj()],
#'    [qwalkr::get_eigschur.spectral()]
#' @export
#' @examples
#' s <- spectral(rbind(c(0.5, 0.3), c(0.3,0.7)))
#'
#' get_eigschur(s, 1, 2) #-> get_eigschur.spectral(...)
#'
get_eigschur <- function(object, ...) UseMethod("get_eigschur")

#' Extract a Schur Cross-Product from a Hermitian Matrix
#'
#' @description
#' Get the Schur product between eigen-projectors based  on the representation of a
#' Hermitian Matrix given by class `spectral`.
#'
#' @param object an instance of class `spectral`.
#' @param id1 index for the first eigenspace according to the ordered (decreasing) spectra.
#' @param id2 index for the second eigenspace according to the ordered (decreasing) spectra. If not provided,
#'    it takes the same value as `id1`.
#' @param ... further arguments passed to or from other methods.
#'
#' @returns The Schur product of the corresponding eigenprojectors, \eqn{E_{id_1} \circ E_{id_2}}.
#'
#' @seealso [qwalkr::spectral()], [qwalkr::get_eigschur()]
#'
#' @export
#'
#' @examples
#' # Spectra is {2, -1} with multiplicities one and two respectively.
#' decomp <- spectral(matrix(c(0,1,1,1,0,1,1,1,0), nrow=3))
#'
#' # Returns the Schur product between the 2-projector and -1-projector.
#' get_eigschur(decomp, id1=2, id2=1)
#'
#' # Returns the Schur square of the 2-projector.
#' get_eigschur(decomp, id1=1, id2=1)
#'
#' # Also returns the Schur square of the 2-projector
#' get_eigschur(decomp, id1=1)
#'
get_eigschur.spectral <- function(object, id1, id2=NULL, ...){

  id2 <- if(is.null(id2)) id1 else id2

  if (out_of_bounds(c(id1, id2), 1, length(object$multiplicity))){
    stop("Index out of bounds! Check the length of the spectra.")
  }
  E_r <- get_eigproj.spectral(object, id1)
  E_s <- get_eigproj.spectral(object, id2)
  return (E_r * E_s)
}

#' Apply a Function to an Operator
#'
#' @param object a representation of the operator.
#' @param ... further arguments passed to or from other methods.
#'
#' @returns The resulting operator from the application of the function.
#' @seealso [qwalkr::act_eigfun.spectral()]
#' @export
#' @examples
#' s <- spectral(rbind(c(0.5, 0.3), c(0.3,0.7)))
#'
#' act_eigfun(s, function(x) x^2) #-> act_eigfun.spectral(...)
#'
act_eigfun <- function(object, ...) UseMethod("act_eigfun")

#' Apply a Function to a Hermitian Matrix
#'
#' @description Apply a function to a Hermitian matrix based on the representation
#' given by class `spectral`.
#'
#' @param object an instance of class `spectral`.
#' @param FUN the function to be applied to the matrix.
#' @param ... further arguments passed on to `FUN`.
#'
#' @returns The matrix resulting from the application of `FUN`.
#'
#'   A Hermitian Matrix admits the spectral decomposition
#'   \deqn{H = \sum_k \lambda_k E_k}
#'   where \eqn{\lambda_k} are its eigenvalues and \eqn{E_k} the
#'   orthogonal projector onto the \eqn{\lambda_k}-eigenspace.
#'
#'   If \eqn{f}=`FUN` is defined on the eigenvalues of `H`, then
#'   `act_eigfun` performs the following calculation
#'
#'   \deqn{f(H) = \sum_k f(\lambda_k) E_k}
#'
#' @seealso [qwalkr::spectral()], [qwalkr::act_eigfun()]
#' @export
#'
#' @examples
#' H <- matrix(c(0,1,1,1,0,1,1,1,0), nrow=3)
#' decomp <- spectral(H)
#'
#' # Calculates H^2.
#' act_eigfun(decomp, FUN = function(x) x^2)
#'
#' # Calculates sin(H).
#' act_eigfun(decomp, FUN = function(x) sin(x))
#'
#' # Calculates H^3.
#' act_eigfun(decomp, FUN = function(x, y) x^y, 3)
#'
act_eigfun.spectral <- function(object, FUN, ...){
  FUN <- match.fun(FUN)

  flam <- FUN(rep(object$eigvals, times=object$multiplicity), ...)

  V <- object$eigvectors

  fH <- V %*% (Conj(t(V)) * flam)
  return (fH)

}

