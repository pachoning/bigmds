% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gower_interpolation_mds.R
\name{gower_interpolation_mds}
\alias{gower_interpolation_mds}
\title{MDS based on Gower interpolation formula}
\usage{
gower_interpolation_mds(x, l, k, dist_fn = stats::dist, ...)
}
\arguments{
\item{x}{A matrix with n individuals (rows) and q variables (columns).}

\item{l}{The largest value which allows classical MDS to be computed efficiently, i.e, the largest value which makes
\code{cmdscale()} be run without any computational issues.}

\item{k}{Number of principal coordinates to be extracted.}

\item{dist_fn}{Distance function to be used for obtaining a MDS configuration.}

\item{...}{Further arguments passed to \code{dist_fn} function.}
}
\value{
Returns a list containing the following elements:
\describe{
\item{points}{A matrix that consists of n individuals (rows) and \code{k} variables (columns) corresponding to the
MDS coordinates.}
\item{eigen}{The first \code{k} eigenvalues.}
}
}
\description{
Performs \emph{Multidimensional Scaling} for big datasets using Gower interpolation formula. This method can
compute a MDS configuration even when the dataset is so large that classical MDS methods (\code{cmdscale}) can not be run
due to computational problems.
}
\details{
\emph{Gower interpolation formula} is the central piece of this algorithm since it allows to add a new set of
points to an existing MDS configuration so that the new one has the same coordinate system.

Given the matrix \code{x} with n individuals (rows) and q variables (columns), a submatrix based on a random sample
of \code{l} individuals is taken and it is used to compute a MDS configuration.

The remaining part of \code{x} is divided into p=(n-\code{l})/\code{l} submatrices. For every submatrix, it is obtained
a MDS configuration by means of \emph{Gower interpolation formula} and the first (random) submatrix. Every MDS
configuration is appended to the existing one so that, at the end of the process, a MDS configuration for \code{x} is
built.
}
\examples{
set.seed(42)
x <- matrix(data = rnorm(4*10000), nrow = 10000) \%*\% diag(c(15, 10, 1, 1))
mds <- gower_interpolation_mds(x = x, l = 200, k = 2, dist_fn = stats::dist)
head(cbind(mds$points, x[, 1:2]))
var(x)
var(mds$points)
}
\references{
Gower, J.C. and D.J, Hand (1995). \emph{Biplots}. Volume 54. CRC Press.

Borg, I. and Groenen, P. (2005). \emph{Modern Multidimensional Scaling: Theory and Applications}. Springer.
}
