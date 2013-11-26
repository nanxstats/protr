\name{extractAAC}
\alias{extractAAC}
\title{Amino Acid Composition Descriptor}
\usage{
  extractAAC(x)
}
\arguments{
  \item{x}{A character vector, as the input protein
  sequence.}
}
\value{
  A length 20 named vector
}
\description{
  Amino Acid Composition Descriptor
}
\details{
  This function calculates the Amino Acid Composition
  descriptor (Dim: 20).
}
\examples{
x = readFASTA(system.file('protseq/P00750.fasta', package = 'protr'))[[1]]
extractAAC(x)
}
\author{
  Nan Xiao <\url{http://www.road2stat.com}>
}
\references{
  M. Bhasin, G. P. S. Raghava. Classification of Nuclear
  Receptors Based on Amino Acid Composition and Dipeptide
  Composition. \emph{Journal of Biological Chemistry},
  2004, 279, 23262.
}
\seealso{
  See \code{\link{extractDC}} and \code{\link{extractTC}}
  for Dipeptide Composition and Tripeptide Composition
  descriptors.
}
\keyword{AAC}
\keyword{Acid}
\keyword{Amino}
\keyword{Composition}
\keyword{extract}
\keyword{extractAAC}

