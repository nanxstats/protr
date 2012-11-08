\name{extractMoreauBroto}
\alias{extractMoreauBroto}
\title{Normalized Moreau-Broto Autocorrelation Descriptor}
\usage{
  extractMoreauBroto(x,
    props = c("ANDN920101", "ARGP820101", "ARGP820102", "ARGP820103", "BEGF750101", "BEGF750102", "BEGF750103", "BHAR880101"),
    nlag = 30L)
}
\arguments{
  \item{x}{A character vector, as the input protein
  sequence.}

  \item{props}{A character vector, specifying the Accession
  Number of the target properties.  8 properties are used
  by default.}

  \item{nlag}{Maximum value of the lag parameter. Default
  is \code{30}.}
}
\value{
  A length \code{nlag} named vector
}
\description{
  Normalized Moreau-Broto Autocorrelation Descriptor
}
\details{
  This function calculates the normalized Moreau-Broto
  autocorrelation descriptor (total \code{length(props) *
  nlag}).
}
\examples{
A06852 = readFASTA(system.file('AAseq/A06852.fasta', package = 'rdpi'))
x = 'MPRLFSYLLGVWLLLSQLPREIPGQSTNDFIKACGRELVRLWVEICGSVSWGRTALSLEEPQLETGPPAETMPSSITKDAEILKMMLEFVPNLPQELKATLSERQPSLRELQQSASKDSNLNFEEFKKIILNRQNEAEDKSLLELKNLGLDKHSRKKRLFRMTLSEKCCQVGCIRKDIARLC'
extractMoreauBroto(x)
}
\author{
  Xiao Nan <\url{http://www.road2stat.com}>
}
\references{
  References needed here!
}
\seealso{
  See \code{\link{extractMoran}} and
  \code{\link{extractGeary}} for Moran autocorrelation
  descriptors and Geary autocorrelation descriptors.
}
\keyword{autocorrelation}
\keyword{Broto}
\keyword{extract}
\keyword{extractMoreauBroto}
\keyword{Moreau}
\keyword{Moreau-Broto}
\keyword{MoreauBroto}
\keyword{normalized}

