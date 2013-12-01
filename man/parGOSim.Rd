\name{parGOSim}
\alias{parGOSim}
\title{Parallellized Protein Sequence Similarity Calculation based on Gene Ontology (GO) Similarity}
\usage{
  parGOSim(x)
}
\arguments{
  \item{x}{A character vector, as the input protein
  sequence.}
}
\value{
  A length \code{lag * p^2} named vector, \code{p} is the
  number of scales (principal components) selected.
}
\description{
  Parallellized Protein Sequence Similarity Calculation
  based on Gene Ontology (GO) Similarity
}
\details{
  This function implemented the parallellized version for
  calculating protein sequence similarity based on Gene
  Ontology (GO) similarity.
}
\note{
  Note
}
\examples{
parGOSim(x)
}
\author{
  Nan Xiao <\url{http://www.road2stat.com}>
}
\seealso{
  See \code{\link{parSeqSim}} for paralleled protein
  similarity calculation based on Smith-Waterman local
  alignment.
}
\keyword{GO}
\keyword{Gene}
\keyword{Ontology}
\keyword{parGOSim}
\keyword{parallel}
\keyword{similarity}

