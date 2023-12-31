% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/linalg.R
\name{get_eigproj.spectral}
\alias{get_eigproj.spectral}
\title{Extract an Eigen-Projector from a Hermitian Matrix}
\usage{
\method{get_eigproj}{spectral}(object, id, ...)
}
\arguments{
\item{object}{an instance of class \code{spectral}.}

\item{id}{index for the desired eigenspace according to the ordered (decreasing) spectra.}

\item{...}{further arguments passed to or from other methods.}
}
\value{
The orthogonal projector of the desired eigenspace.

A Hermitian matrix \code{S} admits the spectral decomposition \eqn{S = \sum_{r}\lambda_r E_r}
such that \eqn{E_r} is the orthogonal projector onto the \eqn{\lambda_r}-eigenspace. If \eqn{V_{id}}
is the matrix associated to the eigenspace, then

\deqn{E_{id} = V_{id}V_{id}^*}
}
\description{
Get the orthogonal projector associated with an eigenspace based on the representation
of a Hermitian Matrix given by class \code{spectral}.
}
\examples{
# Spectra is {2, -1} with multiplicities one and two respectively.
decomp <- spectral(matrix(c(0,1,1,1,0,1,1,1,0), nrow=3))

# Returns the projector associated to the eigenvalue -1.
get_eigproj(decomp, id=2)

# Returns the projector associated to the eigenvalue 2.
get_eigproj(decomp, id=1)

}
\seealso{
\code{\link[=spectral]{spectral()}}, \code{\link[=get_eigproj]{get_eigproj()}}
}
