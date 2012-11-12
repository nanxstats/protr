\name{extractTC}
\alias{extractTC}
\title{Tripeptide Composition Descriptor}
\usage{
  extractTC(x)
}
\arguments{
  \item{x}{A character vector, as the input protein
  sequence.}
}
\value{
  A length 8000 named vector
}
\description{
  Tripeptide Composition Descriptor
}
\details{
  This function calculates the Tripeptide Composition
  descriptor (Dim: 8000).
}
\examples{
x = readFASTA(system.file('protseq/P00750.fasta', package = 'rdpi'))[[1]]
extractTC(x)
}
\author{
  Xiao Nan <\url{http://www.road2stat.com}>
}
\references{
  M. Bhasin, G. P. S. Raghava. Classification of Nuclear
  Receptors Based on Amino Acid Composition and Dipeptide
  Composition. \emph{Journal of Biological Chemistry},
  2004, 279, 23262.
}
\seealso{
  See \code{\link{extractAAC}} and \code{\link{extractDC}}
  for Amino Acid Composition and Dipeptide Composition
  descriptors.
}
\keyword{Composition}
\keyword{extract}
\keyword{extractTC}
\keyword{TC}
\keyword{Tripeptide}

