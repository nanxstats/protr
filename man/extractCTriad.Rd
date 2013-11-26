\name{extractCTriad}
\alias{extractCTriad}
\title{Conjoint Triad Descriptor}
\usage{
  extractCTriad(x)
}
\arguments{
  \item{x}{A character vector, as the input protein
  sequence.}
}
\value{
  A length 343 named vector
}
\description{
  Conjoint Triad Descriptor
}
\details{
  This function calculates the Conjoint Triad descriptor
  (Dim: 343).
}
\examples{
x = readFASTA(system.file('protseq/P00750.fasta', package = 'protr'))[[1]]
extractCTriad(x)
}
\author{
  Nan Xiao <\url{http://www.road2stat.com}>
}
\references{
  J.W. Shen, J. Zhang, X.M. Luo, W.L. Zhu, K.Q. Yu, K.X.
  Chen, Y.X. Li, H.L. Jiang. Predicting Protein-protein
  Interactions Based Only on Sequences Information.
  \emph{Proceedings of the National Academy of Sciences}.
  007, 104, 4337-4341.
}
\keyword{CTriad}
\keyword{Conjoint}
\keyword{Triad}
\keyword{extract}
\keyword{extractCTriad}

