

#' Create a Continuous-time Quantum Walk
#'
#' @description Function to create a continuous-time quantum walk (ctqwalk) object.
#'
#' @param hamiltonian  a Hermitian Matrix representing the Hamiltonian of the system.
#'
#' @return An object of class "ctqwalk" containing the continuous-time quantum walk
#'   object.
#' @export
#'
#' @examples
#' # Quantum walk on the Path graph with three vertices, the hamiltonian is the adjacency matrix.
#' P3 <- matrix(c(0,1,0,1,0,1,0,1,0), nrow=3)
#'
#' walk <- ctqwalk(P3)
#'
#'
ctqwalk <- function(hamiltonian){
  if(!isSymmetric(hamiltonian)){
    warning("The Hamiltonian must be a Hermitian matrix!")
    return(NULL)
  }
  spectral <- eigen(hamiltonian, symmetric = TRUE)

  output <- list(hamiltonian = hamiltonian,
                 eigvals=spectral$values,
                 eigid = label_unique(spectral$values),
                 eigvectors=spectral$vectors)


  class(output) <- c("ctqwalk", "qwalk")

  return(output)
}



