#' Get the Average Mixing Matrix of a Quantum Walk
#'
#' @param object a representation of the quantum walk.
#' @param ... further arguments passed to or from other methods.
#'
#' @returns The average mixing matrix.
#' @export
#'
#' @seealso [qwalkr::mixing_matrix], [qwalkr::avg_matrix.ctqwalk]
avg_matrix <- function(object, ...) UseMethod("avg_matrix")


#' Get the Average Mixing Matrix of a Continuous-Time Quantum Walk
#'
#' @param object a representation of the quantum walk.
#' @param ... further arguments passed to or from other methods.
#'
#' @returns The average mixing matrix of the continuous-time quantum walk.
#'
#'    Let \eqn{M(t)} be the mixing matrix of the quantum walk, then the
#'    average mixing matrix is defined as
#'
#'    \deqn{\widehat{M} := \lim_{T \to \infty} \frac{1}{T}\int_{0}^T M(t)\textrm{d}t}
#'
#'    and encodes the long-term average behavior of the walk. It is possible to prove that
#'
#'    \deqn{\widehat{M} = \sum_r E_r \circ E_r}
#'
#'    in which \eqn{E_r} is the orthogonal projector onto the \eqn{\lambda_r}-eigenspace.
#'
#'    `avg_matrix` returns the matrix from the above calculation.
#'
#' @export
#'
#' @seealso [qwalkr::ctqwalk], [qwalkr::avg_matrix]
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



