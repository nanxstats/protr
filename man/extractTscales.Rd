\name{extractTscales}
\alias{extractTscales}
\title{T-scales Descriptor}
\usage{
  extractTscales(x)
}
\arguments{
  \item{x}{A character vector, as the input protein
  sequence.}
}
\value{
  A length \code{20 + n * lambda} named vector, \code{n} is
  the number of properties selected.
}
\description{
  T-scales Descriptor
}
\details{
  This function calculates the T-scales descriptor (Dim:
  \code{20 + (n * lambda)}, \code{n} is the number of
  properties selected, default is 80).
}
\note{
  Note
}
\examples{
x = readFASTA(system.file('protseq/P00750.fasta', package = 'protr'))[[1]]
extractTscales(x)
}
\author{
  Xiao Nan <\url{http://www.road2stat.com}>
}
\references{
  TBD
}
\seealso{
  See \code{\link{extractZscales}},
  \code{\link{extractSTscales}} and
  \code{\link{extractVHSE}} for Z-scales, ST-scales and
  VHSE descriptors. See \code{\link{extractScales}} for
  generalized scales-based descriptors.
}
\keyword{extract}
\keyword{extractTscales}
\keyword{PCM}
\keyword{T-scales}

