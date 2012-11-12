\name{extractDC}
\alias{extractDC}
\title{Dipeptide Composition Descriptor}
\usage{
  extractDC(x)
}
\arguments{
  \item{x}{A character vector, as the input protein
  sequence.}
}
\value{
  A length 400 named vector
}
\description{
  Dipeptide Composition Descriptor
}
\details{
  This function calculates the Dipeptide Composition
  descriptor (Dim: 400).
}
\examples{
x = readFASTA(system.file('protseq/P00750.fasta', package = 'rdpi'))[[1]]
extractDC(x)
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
  See \code{\link{extractAAC}} and \code{\link{extractTC}}
  for Amino Acid Composition and Tripeptide Composition
  descriptors.
}
\keyword{Composition}
\keyword{DC}
\keyword{Dipeptide}
\keyword{extract}
\keyword{extractDC}

