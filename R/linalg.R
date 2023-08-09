#' Spectral Decomposition of a Matrix
#'
#' @description `spectral()` is a wrapper around [base::eigen()] which can handle repeated eigenvalues.
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
#' @param id index of the eigenspace according to the ordered (decreasing) spectra.
#' @param ... further arguments passed to or from other methods.
#' @returns A matrix whose columns form the orthonormal eigenbasis.
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
#' @param id index of the eigenspace according to the ordered (decreasing) spectra.
#' @param ... further arguments passed to or from other methods.
#'
#' @returns The orthogonal projector of the desired eigenspace. Note
#'    that its rank is equal to the multiplicity of the eigenvalue.
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
  return (A %*% t(A))
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
#' @param id1 index of the first eigenspace according to the ordered (decreasing) spectra.
#' @param id2 index of the second eigenspace according to the ordered (decreasing) spectra. Defaults to
#'   the same index as id1.
#' @param ... further arguments passed to or from other methods.
#'
#' @returns The Schur product of the corresponding eigenprojectors.
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
#' # Also returns the Schur square of the 2-projector (Default of id2 is id1)
#' extractSCHUR(decomp, id1=1)
#'
extractSCHUR.spectral <- function(object, id1, id2=id1, ...){
  if (out_of_bounds(c(id1, id2), 1, length(object$multiplicity))){
    warning("Index out of bounds! Check the length of the spectra.")
    return(NULL)
  }
  E_r <- extractPROJ.spectral(object, id1)
  E_s <- extractPROJ.spectral(object, id2)
  return (E_r * E_s)
}



