\name{extractCTDD}
\alias{extractCTDD}
\title{CTD Descriptors - Distribution}
\usage{
  extractCTDD(x)
}
\arguments{
  \item{x}{A character vector, as the input protein
  sequence.}
}
\value{
  A length 105 named vector
}
\description{
  CTD Descriptors - Distribution
}
\details{
  This function calculates the Distribution descriptor of
  the CTD descriptors (Dim: 105).
}
\examples{
x = readFASTA(system.file('protseq/P00750.fasta', package = 'rdpi'))[[1]]
extractCTDD(x)
}
\author{
  Xiao Nan <\url{http://www.road2stat.com}>
}
\references{
  Inna Dubchak, Ilya Muchink, Stephen R. Holbrook and
  Sung-Hou Kim. Prediction of protein folding class using
  global description of amino acid sequence.
  \emph{Proceedings of the National Academy of Sciences}.
  USA, 1995, 92, 8700-8704.

  Inna Dubchak, Ilya Muchink, Christopher Mayor, Igor
  Dralyuk and Sung-Hou Kim. Recognition of a Protein Fold
  in the Context of the SCOP classification.
  \emph{Proteins: Structure, Function and Genetics}, 1999,
  35, 401-407.
}
\seealso{
  See \code{\link{extractCTDC}} and
  \code{\link{extractCTDT}} for Composition and Transition
  of the CTD descriptors.
}
\keyword{Composition}
\keyword{CTD}
\keyword{CTDD}
\keyword{extract}
\keyword{extractCTDD}

