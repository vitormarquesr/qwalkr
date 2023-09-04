#' Create a Continuous-time Quantum Walk
#'
#' @description `ctqwalk()` creates a quantum walk object from a hamiltonian.
#'
#' @param hamiltonian  a Hermitian Matrix representing the Hamiltonian of the system.
#' @param ... further arguments passed on to [qwalkr::spectral()]
#'
#' @returns A list with the walk related objects, i.e the hamiltonian and its spectral
#'   decomposition (See [qwalkr::spectral()] for further details)
#'
#' @seealso [qwalkr::spectral()], [qwalkr::unitary_matrix.ctqwalk],
#' [qwalkr::mixing_matrix.ctqwalk], [qwalkr::avg_matrix.ctqwalk],
#' [qwalkr::gavg_matrix.ctqwalk]
#' @export
#'
#' @examples
#' # Creates a walk from the adjacency matrix of the graph P3.
#' ctqwalk(matrix(c(0,1,0,1,0,1,0,1,0), nrow=3))
#'
#'
ctqwalk <- function(hamiltonian, ...){
  if(!isSymmetric(hamiltonian)){
    stop("Hamiltonian must be a Hermitian matrix!")
  }

  output <- c(list(hamiltonian=hamiltonian), spectral(hamiltonian, ...))

  class(output) <- c("ctqwalk", "spectral")

  return (output)
}

