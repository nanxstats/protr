\name{extractFASGAI}
\alias{extractFASGAI}
\title{Factor Analysis Scales of Generalized AA Information (FASGAI) Descriptor}
\usage{
  extractFASGAI(x)
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
  Factor Analysis Scales of Generalized AA Information
  (FASGAI) Descriptor
}
\details{
  This function calculates the Factor Analysis Scales of
  Generalized AA Information (FASGAI) descriptors. (Dim:
  \code{20 + (n * lambda)}, \code{n} is the number of
  properties selected, default is 80).
}
\note{
  Note
}
\examples{
x = readFASTA(system.file('protseq/P00750.fasta', package = 'protr'))[[1]]
extractFASGAI(x)
}
\author{
  Xiao Nan <\url{http://www.road2stat.com}>
}
\references{
  TBD
}
\keyword{extract}
\keyword{extractFASGAI}
\keyword{FASGAI}
\keyword{PCM}

