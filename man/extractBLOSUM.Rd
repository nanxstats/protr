\name{extractBLOSUM}
\alias{extractBLOSUM}
\title{BLOSUM Descriptor}
\usage{
  extractBLOSUM(x)
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
  BLOSUM Descriptor
}
\details{
  This function calculates the BLOSUM descriptor (Dim:
  \code{20 + (n * lambda)}, \code{n} is the number of
  properties selected, default is 80).
}
\note{
  Note
}
\examples{
x = readFASTA(system.file('protseq/P00750.fasta', package = 'protr'))[[1]]
# extractBLOSUM(x)
}
\author{
  Nan Xiao <\url{http://www.road2stat.com}>
}
\references{
  TBD
}
\keyword{BLOSUM}
\keyword{PCM}
\keyword{extract}
\keyword{extractBLOSUM}

