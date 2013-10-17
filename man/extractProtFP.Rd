\name{extractProtFP}
\alias{extractProtFP}
\title{Protein Fingerprint (ProtFP)}
\usage{
  extractProtFP(x)
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
  Protein Fingerprint (ProtFP)
}
\details{
  This function calculates the Protein Fingerprint (ProtFP)
  (Dim: \code{20 + (n * lambda)}, \code{n} is the number of
  properties selected, default is 80).
}
\note{
  Note
}
\examples{
x = readFASTA(system.file('protseq/P00750.fasta', package = 'protr'))[[1]]
extractProtFP(x)
}
\author{
  Xiao Nan <\url{http://www.road2stat.com}>
}
\references{
  TBD
}
\keyword{extract}
\keyword{extractProtFP}
\keyword{fingerprint}
\keyword{PCM}
\keyword{ProtFP}

