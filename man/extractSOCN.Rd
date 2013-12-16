\name{extractSOCN}
\alias{extractSOCN}
\title{Sequence-Order-Coupling Numbers}
\usage{
extractSOCN(x, nlag = 30)
}
\arguments{
  \item{x}{A character vector, as the input protein
  sequence.}

  \item{nlag}{The maximum lag, defualt is 30.}
}
\value{
A length \code{nlag * 2} named vector
}
\description{
Sequence-Order-Coupling Numbers
}
\details{
This function calculates the Sequence-Order-Coupling
Numbers (Dim: \code{nlag * 2}, default is 60).
}
\examples{
x = readFASTA(system.file('protseq/P00750.fasta', package = 'protr'))[[1]]
extractSOCN(x)
}
\author{
Nan Xiao <\url{http://www.road2stat.com}>
}
\references{
Kuo-Chen Chou. Prediction of Protein Subcellar Locations by
Incorporating Quasi-Sequence-Order Effect.
\emph{Biochemical and Biophysical Research Communications},
2000, 278, 477-483.

Kuo-Chen Chou and Yu-Dong Cai. Prediction of Protein
Sucellular Locations by GO-FunD-PseAA Predictor.
\emph{Biochemical and Biophysical Research Communications},
2004, 320, 1236-1239.

Gisbert Schneider and Paul Wrede. The Rational Design of
Amino Acid Sequences by Artifical Neural Networks and
Simulated Molecular Evolution: Do Novo Design of an
Idealized Leader Cleavge Site. \emph{Biophys Journal},
1994, 66, 335-344.
}
\seealso{
See \code{\link{extractQSO}} for quasi-sequence-order
descriptors.
}
\keyword{Coupling}
\keyword{Number}
\keyword{Order}
\keyword{SOCN}
\keyword{Sequence}
\keyword{extract}
\keyword{extractSOCN}

