\name{extractVHSE}
\alias{extractVHSE}
\title{Vectors of Hydrophobic, Steric, and Electronic (VHSE) Descriptor}
\usage{
  extractVHSE(x)
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
  Vectors of Hydrophobic, Steric, and Electronic (VHSE)
  Descriptor
}
\details{
  This function calculates the Vectors of Hydrophobic,
  Steric, and Electronic (VHSE) descriptor. (Dim: \code{20
  + (n * lambda)}, \code{n} is the number of properties
  selected, default is 80).
}
\note{
  Note
}
\examples{
x = readFASTA(system.file('protseq/P00750.fasta', package = 'protr'))[[1]]
extractVHSE(x)
}
\author{
  Xiao Nan <\url{http://www.road2stat.com}>
}
\references{
  TBD
}
\seealso{
  See \code{\link{extractZscales}},
  \code{\link{extractTscales}} and
  \code{\link{extractSTscales}} for Z-scales, T-scales and
  ST-scales descriptors. See \code{\link{extractScales}}
  for generalized scales-based descriptors.
}
\keyword{extract}
\keyword{extractVHSE}
\keyword{PCM}
\keyword{VHSE}

