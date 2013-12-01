\name{parSeqSim}
\alias{parSeqSim}
\title{Parallellized Protein Sequence Similarity Calculation based on Sequence Alignment}
\usage{
  parSeqSim(x)
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
  based on Sequence Alignment
}
\details{
  This function implemented the parallellized version for
  calculating protein sequence similarity based on sequence
  alignment.
}
\note{
  Note
}
\examples{
parSeqSim(x)
}
\author{
  Nan Xiao <\url{http://www.road2stat.com}>
}
\seealso{
  See \code{\link{parGOSim}} for paralleled protein
  similarity calculation based on Gene Ontology (GO)
  similarity.
}
\keyword{alignment}
\keyword{local}
\keyword{parSeqSim}
\keyword{parallel}
\keyword{sequence}
\keyword{similarity}
\keyword{smith-waterman}

