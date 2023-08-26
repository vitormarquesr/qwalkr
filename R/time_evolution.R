
#' Get the Unitary Time Evolution Operator of a Quantum Walk
#'
#' @param object a representation of the quantum walk.
#' @param ... further arguments passed to or from other methods.
#'
#' @returns The unitary time evolution operator.
#' @export
#'
#' @seealso [qwalkr::unitary_matrix.ctqwalk]
#'
unitary_matrix <- function(object, ...) UseMethod("unitary_matrix")


#' Get the Unitary Time Evolution Operator of a Continuous-Time Quantum Walk
#'
#' @param object an instance of class `ctqwalk`.
#' @param t it will be returned the evolution operator at time `t`.
#' @param ... further arguments passed to or from other methods.
#'
#' @returns The unitary time evolution operator of the continuous-time quantum walk
#'   at time `t`.
#'
#'   If \eqn{|\psi(t) \rangle} is the quantum state of the system at time \eqn{t}, and
#'   \eqn{H} the Hamiltonian operator, then the evolution is governed by
#'   the Schrodinger equation
#'
#'   \deqn{H|\psi(t) \rangle = i \frac{\partial}{\partial t}|\psi(t) \rangle}
#'
#'    and if \eqn{H} is time-independent its solution is given by
#'
#'    \deqn{|\psi(t) \rangle = U(t)|\psi(0) \rangle = e^{-iHt}|\psi(0) \rangle}
#'
#'    The evolution operator is the result of the complex matrix exponential
#'    and it can be calculated as
#'
#'    \deqn{U(t) = e^{-iHt} = \sum_r e^{-i t \lambda_r}E_r}
#'
#'    if \eqn{H = \sum_r \lambda_r E_r}.
#'
#'    `unitary_matrix` returns the result of the above calculation for a
#'    provided `t`.
#'
#' @export
#' @seealso [qwalkr::ctqwalk], [qwalkr::unitary_matrix], [qwalkr::act_eigfun]
#'
#' @examples
#' walk <- ctqwalk(matrix(c(0,1,0,1,0,1,0,1,0), nrow=3))
#'
#' # Returns the operator at time t = 2*pi, U(2pi)
#' unitary_matrix(walk, t = 2*pi)
#'
unitary_matrix.ctqwalk <- function(object, t, ...){
  U <- act_eigfun(object, function(x) exp(-1i*x*t))
  return (U)
}





