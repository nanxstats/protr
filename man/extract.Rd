\name{extract}
\alias{extract}
\title{Feature Extraction for Protein Sequences}
\usage{
  extract(x, type, ...)
}
\arguments{
  \item{x}{A character vector, as the input protein
  sequence.  Length \code{n}.}

  \item{type}{a vector of size n * 1 describing the
  chunklets: -1 in the i'th place says that point i doesn't
  belong to any chunklet; integer j in place i says that
  point i belongs to chunklet j.  The chunklets indexes
  should be 1:(number of chunklets).}

  \item{\dots}{Additional options used by specific
  descriptor as described below.}
}
\value{
  The result data frame

  The three returned argument are just different forms of
  the same output. If one is interested in a Mahalanobis
  metric over the original data space, the first argument
  is all she/he needs. If a transformation into another
  space (where one can use the Euclidean metric) is
  preferred, the second returned argument is sufficient.
  Using A and B is equivalent in the following sense:

  \describe{ \item{\code{'ReliefFequalK'}}{ReliefF
  algorithm where k nearest instances have equal weight.
  ReliefF algorithm where k nearest instances have equal
  weight.} \item{\code{'ReliefFequalK'}}{ReliefF algorithm
  where k nearest instances have equal weight.}
  \item{\code{'ReliefFequalK'}}{ReliefF algorithm where k
  nearest instances have equal weight.} }
}
\description{
  Feature Extraction for Protein Sequences
}
\details{
  This function takes a protein sequence and descriptor
  type as arguments and returns the result as a data frame.
}
\note{
  Note that any different sets of instances (chunklets),
  e.g. {1, 3, 7} and {4, 6}, might belong to the same class
  and might belong to different classes.
}
\examples{
require(seqinr)
AAseq = read.fasta(system.file('AAseq/AAseq.fasta', package = 'rdpi'))
extract(AAseq, 'AAC')
}
\author{
  Xiao Nan <\url{http://www.road2stat.com}>
}
\references{
  Aharon Bar-Hillel, Tomer Hertz, Noam Shental, and Daphna
  Weinshall (2003). Learning Distance Functions using
  Equivalence Relations. \emph{Proceedings of 20th
  International Conference on Machine Learning (ICML2003)}.
}
\seealso{
  See \code{\link{extractAAC}} for Amino Acid Composition,
  etc.
}
\keyword{descriptors}
\keyword{extract}
\keyword{extraction}
\keyword{feature}
\keyword{rdpi}

