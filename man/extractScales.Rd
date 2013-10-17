\name{extractScales}
\alias{extractScales}
\title{Generalized Scales-Based Descriptors}
\usage{
  extractScales(x, propmat, pc, lag, scale = TRUE)
}
\arguments{
  \item{x}{A character vector, as the input protein
  sequence.}

  \item{propmat}{A matrix containing the properties for the
  amino acids.  Each row represent one amino acid type,
  each column represents one property.  Note that the
  one-letter row names must be provided for we need them to
  seek the properties for each AA type.}

  \item{pc}{Integer. Use the first pc principle components
  as the scales.  Must be no greater than the number of AA
  properties provided.}

  \item{lag}{The lag parameter. Must be less than the amino
  acids.}

  \item{scale}{Logical. Should we auto-scale the property
  matrix (\code{propmat}) before PCA?}
}
\value{
  A length \code{lag * p^2} named vector, \code{p} is the
  number of scales (principle components) selected.
}
\description{
  Generalized Scales-Based Descriptors
}
\details{
  This function calculates the generalized scales-based
  descriptors. This function implements the core
  computation needed for Z-scales T-scales, ST-scales and
  VHSE descriptors in the protr package.
}
\note{
  Note
}
\examples{
x = readFASTA(system.file('protseq/P00750.fasta', package = 'protr'))[[1]]
extractScales(x)
}
\author{
  Xiao Nan <\url{http://www.road2stat.com}>
}
\seealso{
  See \code{\link{extractZscales}},
  \code{\link{extractTscales}},
  \code{\link{extractSTscales}} and
  \code{\link{extractVHSE}} for Z-scales, T-scales,
  ST-scales and VHSE descriptors.
}
\keyword{extract}
\keyword{extractScales}
\keyword{PCM}
\keyword{scales}

