#' The Average Mixing Matrix of a Quantum Walk
#'
#' @param object a representation of the quantum walk.
#' @param ... further arguments passed to or from other methods.
#'
#' @returns The average mixing matrix.
#' @export
#'
#' @seealso [qwalkr::mixing_matrix()], [qwalkr::gavg_matrix()], [qwalkr::avg_matrix.ctqwalk()]
#' @examples
#' w <- ctqwalk(matrix(c(0,1,0,1,0,1,0,1,0), nrow=3))
#'
#' avg_matrix(w) #-> avg_matrix.ctqwalk(...)
#'
avg_matrix <- function(object, ...) UseMethod("avg_matrix")


#' The Average Mixing Matrix of a Continuous-Time Quantum Walk
#'
#' @details Let \eqn{M(t)} be the mixing matrix of the quantum walk, then the average mixing matrix is defined as
#'
#'    \deqn{\widehat{M} := \lim_{T \to \infty} \frac{1}{T}\int_{0}^T M(t)\textrm{d}t}
#'
#'    and encodes the long-term average behavior of the walk. Given the Hamiltonian
#'    \eqn{H = \sum_r \lambda_r E_r}, it is possible to prove that
#'
#'    \deqn{\widehat{M} = \sum_r E_r \circ E_r}
#'
#' @param object a representation of the quantum walk.
#' @param ... further arguments passed to or from other methods.
#'
#' @returns `avg_matrix()` returns the average mixing matrix
#'   as a square matrix of the same order as the walk.
#'
#' @export
#'
#' @seealso [qwalkr::ctqwalk()], [qwalkr::avg_matrix()]
#' @examples
#' walk <- ctqwalk(matrix(c(0,1,0,1,0,1,0,1,0), nrow=3))
#'
#' # Return the average mixing matrix
#' avg_matrix(walk)
avg_matrix.ctqwalk <- function(object, ...){
  ids <- seq_along(object$multiplicity)

  E2 <- lapply(ids, get_eigschur, object=object)

  M <- Reduce(`+`, E2)

  return (M)
}

#' The Generalized Average Mixing Matrix of a Quantum Walk
#'
#' @param object a representation of the quantum walk.
#' @param ... further arguments passed to or from other methods.
#'
#' @returns The generalized average mixing matrix.
#' @export
#'
#' @seealso [qwalkr::mixing_matrix()], [qwalkr::avg_matrix()], [qwalkr::gavg_matrix.ctqwalk()]
#'
#' @examples
#' w <- ctqwalk(matrix(c(0,1,0,1,0,1,0,1,0), nrow=3))
#'
#' gavg_matrix(w, rnorm(100)) #-> gavg_matrix.ctqwalk(...)
#'
gavg_matrix <- function(object, ...) UseMethod("gavg_matrix")


#' The Generalized Average Mixing Matrix of a Continuous-Time Quantum Walk
#'
#' @details Let \eqn{M(t)} be the mixing matrix of the quantum walk and \eqn{R} a random variable
#'   with associated probability density function \eqn{f_R(t)}. Then the generalized average mixing
#'   matrix under \eqn{R} is defined as
#'
#'    \deqn{\widehat{M}_R := \mathbb{E}[M(R)] = \int_{-\infty}^{\infty} M(t)f_R(t)\textrm{d}t}
#'
#' @param object a representation of the quantum walk.
#' @param R samples from the random variable \eqn{R} (For performance, it is recommended at most 10000 samples).
#' @param ... further arguments passed to or from other methods.
#'
#' @returns `gavg_matrix()` returns the generalized average mixing matrix
#'   as a square matrix of the same order as the walk.
#'
#' @export
#'
#' @seealso [qwalkr::ctqwalk()], [qwalkr::gavg_matrix()]
#' @examples
#' walk <- ctqwalk(matrix(c(0,1,0,1,0,1,0,1,0), nrow=3))
#'
#' # Return the average mixing matrix under a Standard Gaussian distribution
#' gavg_matrix(walk, rnorm(1000))
gavg_matrix.ctqwalk <- function(object, R, ...){

  M_t <- lapply(R, mixing_matrix, object=object)

  M_R <- Reduce(`+`, M_t)/length(R)

  return (M_R)

}





