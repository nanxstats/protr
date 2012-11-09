\name{extractSOCN}
\alias{extractSOCN}
\title{Sequence-Order-Coupling Numbers}
\usage{
  extractSOCN(x, dist = c("Schneider-Wrede", "Grantham"),
    nlag = 60)
}
\arguments{
  \item{x}{A character vector, as the input protein
  sequence.}

  \item{dist}{A character, specify the distance matrix used
  between the 20 amini acids. One of 'Schneider-Wrede' or
  'Grantham'}

  \item{nlag}{The maximum lag, defualt is 60.}
}
\value{
  A length \code{nlag} named vector
}
\description{
  Sequence-Order-Coupling Numbers
}
\details{
  This function calculates the Sequence-Order-Coupling
  Numbers (Dim: \code{nlag}).
}
\examples{
A06852 = readFASTA(system.file('AAseq/A06852.fasta', package = 'rdpi'))
x = 'MPRLFSYLLGVWLLLSQLPREIPGQSTNDFIKACGRELVRLWVEICGSVSWGRTALSLEEPQLETGPPAETMPSSITKDAEILKMMLEFVPNLPQELKATLSERQPSLRELQQSASKDSNLNFEEFKKIILNRQNEAEDKSLLELKNLGLDKHSRKKRLFRMTLSEKCCQVGCIRKDIARLC'
extractSOCN(x, dist = 'Schneider-Wrede')
}
\author{
  Xiao Nan <\url{http://www.road2stat.com}>
}
\references{
  Kuo-Chen Chou. Prediction of Protein Subcellar Locations
  by Incorporating Quasi-Sequence-Order Effect.
  \emph{Biochemical and Biophysical Research
  Communications}, 2000, 278, 477-483.

  Gisbert Schneider and Paul Wrede. The Rational Design of
  Amino Acid Sequences by Artificial Neural Networks and
  Simulated Molecular Evolution: Do Novo Design of an
  Idealized Leader Cleavage Site. \emph{Biophys Journal},
  1994, 66, 335-344.

  Grantham, R. Amino acid difference formula to help
  explain protein evolution. \emph{Science}, 1974, 185,
  862-864.
}
\seealso{
  See \code{\link{extractQSO}} for quasi-sequence-order
  descriptors.
}
\keyword{Coupling}
\keyword{extract}
\keyword{extractSOCN}
\keyword{Number}
\keyword{Order}
\keyword{Sequence}
\keyword{SOCN}

