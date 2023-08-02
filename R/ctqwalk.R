

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
#' @seealso [qwalkr::spectral()]
#' @export
#'
#' @examples
#'
#' ctqwalk(matrix(c(0,1,0,1,0,1,0,1,0), nrow=3))
#'
#'
ctqwalk <- function(hamiltonian, ...){
  if(!isSymmetric(hamiltonian)){
    warning("The Hamiltonian must be a Hermitian matrix!")
    return(NULL)
  }

  output <- c(list(hamiltonian=hamiltonian), spectral(hamiltonian, symmetric=TRUE, ...))

  class(output) <- c("ctqwalk", "spectral")

  return (output)
}

