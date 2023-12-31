% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/linalg.R
\name{get_eigschur.spectral}
\alias{get_eigschur.spectral}
\title{Extract a Schur Cross-Product from a Hermitian Matrix}
\usage{
\method{get_eigschur}{spectral}(object, id1, id2 = NULL, ...)
}
\arguments{
\item{object}{an instance of class \code{spectral}.}

\item{id1}{index for the first eigenspace according to the ordered (decreasing) spectra.}

\item{id2}{index for the second eigenspace according to the ordered (decreasing) spectra. If not provided,
it takes the same value as \code{id1}.}

\item{...}{further arguments passed to or from other methods.}
}
\value{
The Schur product of the corresponding eigenprojectors, \eqn{E_{id_1} \circ E_{id_2}}.
}
\description{
Get the Schur product between eigen-projectors based  on the representation of a
Hermitian Matrix given by class \code{spectral}.
}
\examples{
# Spectra is {2, -1} with multiplicities one and two respectively.
decomp <- spectral(matrix(c(0,1,1,1,0,1,1,1,0), nrow=3))

# Returns the Schur product between the 2-projector and -1-projector.
get_eigschur(decomp, id1=2, id2=1)

# Returns the Schur square of the 2-projector.
get_eigschur(decomp, id1=1, id2=1)

# Also returns the Schur square of the 2-projector
get_eigschur(decomp, id1=1)

}
\seealso{
\code{\link[=spectral]{spectral()}}, \code{\link[=get_eigschur]{get_eigschur()}}
}
