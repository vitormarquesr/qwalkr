
#' The Unitary Time Evolution Operator of a Quantum Walk
#'
#' @param object a representation of the quantum walk.
#' @param ... further arguments passed to or from other methods.
#'
#' @returns The unitary time evolution operator.
#' @export
#'
#' @seealso [qwalkr::mixing_matrix()], [qwalkr::unitary_matrix.ctqwalk()]
#'
#' @examples
#' w <- ctqwalk(matrix(c(0,1,0,1,0,1,0,1,0), nrow=3))
#'
#' unitary_matrix(w, t = 2*pi) #-> unitary_matrix.ctqwalk(...)
#'
unitary_matrix <- function(object, ...) UseMethod("unitary_matrix")


#' The Unitary Time Evolution Operator of a Continuous-Time Quantum Walk
#'
#' @details If \eqn{|\psi(t) \rangle} is the quantum state of the system at time \eqn{t}, and
#'   \eqn{H} the Hamiltonian operator, then the evolution is governed by
#'   the Schrodinger equation
#'
#'   \deqn{\frac{\partial}{\partial t}|\psi(t) \rangle = iH|\psi(t) \rangle}
#'
#'    and if \eqn{H} is time-independent its solution is given by
#'
#'    \deqn{|\psi(t) \rangle = U(t)|\psi(0) \rangle = e^{iHt}|\psi(0) \rangle}
#'
#'    The evolution operator is the result of the complex matrix exponential
#'    and it can be calculated as
#'
#'    \deqn{U(t) = e^{iHt} = \sum_r e^{i t \lambda_r}E_r}
#'
#'    in which \eqn{H = \sum_r \lambda_r E_r}.
#'
#' @param object an instance of class `ctqwalk`.
#' @param t it will be returned the evolution operator at time `t`.
#' @param ... further arguments passed to or from other methods.
#'
#' @returns `unitary_matrix()` returns the unitary time evolution operator of the
#'   CTQW evaluated at time `t`.
#'
#' @export
#' @seealso [qwalkr::ctqwalk()], [qwalkr::unitary_matrix()],
#' [qwalkr::act_eigfun()]
#'
#' @examples
#' walk <- ctqwalk(matrix(c(0,1,0,1,0,1,0,1,0), nrow=3))
#'
#' # Returns the operator at time t = 2*pi, U(2pi)
#' unitary_matrix(walk, t = 2*pi)
#'
unitary_matrix.ctqwalk <- function(object, t, ...){
  Ut <- act_eigfun(object, function(x) exp(1i*x*t))
  return (Ut)
}

#' The Mixing Matrix of a Quantum Walk
#'
#' @param object a representation of the quantum walk.
#' @param ... further arguments passed to or from other methods.
#'
#' @returns The mixing matrix of the quantum walk.
#' @seealso [qwalkr::unitary_matrix()], [qwalkr::avg_matrix()], [qwalkr::gavg_matrix()],
#' [qwalkr::mixing_matrix.ctqwalk()]
#' @export
#'
#' @examples
#' w <- ctqwalk(matrix(c(0,1,0,1,0,1,0,1,0), nrow=3))
#'
#' mixing_matrix(w, t = 2*pi) #-> mixing_matrix.ctqwalk(...)
#'
mixing_matrix <- function(object, ...) UseMethod("mixing_matrix")


#' The Mixing Matrix of a Continuous-Time Quantum Walk
#'
#' @details Let \eqn{U(t)} be the time evolution operator of the quantum walk at
#'    time \eqn{t}, then the mixing matrix is given by
#'
#'    \deqn{M(t) = U(t) \circ \overline{U(t)}}
#'
#'    \eqn{M(t)} is a doubly stochastic real symmetric matrix, which encodes the
#'    probability density of the quantum system at time \eqn{t}.
#'
#'    More precisely, the \eqn{(M(t))_{ab}} entry gives us the probability
#'    of measuring the standard basis state \eqn{|b \rangle} at time \eqn{t}, given that
#'    the quantum walk started at \eqn{|a \rangle}.
#'
#' @param object an instance of class `ctqwalk`.
#' @param t it will be returned the mixing matrix at time `t`.
#' @param ... further arguments passed to or from other methods.
#'
#' @returns `mixing_matrix()` returns the  mixing matrix of the CTQW
#'   evaluated at time `t`.
#'
#' @export
#' @seealso [qwalkr::ctqwalk()], [qwalkr::mixing_matrix()]
#'
#' @examples
#' walk <- ctqwalk(matrix(c(0,1,0,1,0,1,0,1,0), nrow=3))
#'
#' # Returns the mixing matrix at time t = 2*pi, M(2pi)
#' mixing_matrix(walk, t = 2*pi)
#'
mixing_matrix.ctqwalk <- function(object, t, ...){
  Mt <- Mod(unitary_matrix(object, t))^2
  return (Mt)
}



