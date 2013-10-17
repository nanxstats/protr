\name{extractWHIM}
\alias{extractWHIM}
\title{MS-WHIM (Molecular Surfaces-Weighted Holistic Invariant Molecular) Descriptor}
\usage{
  extractWHIM(x)
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
  MS-WHIM (Molecular Surfaces-Weighted Holistic Invariant
  Molecular) Descriptor
}
\details{
  This function calculates the MS-WHIM (Molecular
  Surfaces-Weighted Holistic Invariant Molecular)
  descriptors. (Dim: \code{20 + (n * lambda)}, \code{n} is
  the number of properties selected, default is 80).
}
\note{
  Note
}
\examples{
x = readFASTA(system.file('protseq/P00750.fasta', package = 'protr'))[[1]]
extractWHIM(x)
}
\author{
  Xiao Nan <\url{http://www.road2stat.com}>
}
\references{
  TBD
}
\keyword{extract}
\keyword{extractWHIM}
\keyword{PCM}
\keyword{WHIM}

