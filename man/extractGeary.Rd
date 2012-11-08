\name{extractGeary}
\alias{extractGeary}
\title{Geary Autocorrelation Descriptor}
\usage{
  extractGeary(x,
    props = c("CIDH920105", "BHAR880101", "CHAM820101", "CHAM820102", "CHOC760101", "BIGC670101", "CHAM810101", "DAYM780201"),
    nlag = 30L)
}
\arguments{
  \item{x}{A character vector, as the input protein
  sequence.}

  \item{props}{A character vector, specifying the Accession
  Number of the target properties.  8 properties are used
  by default, as listed below: \describe{ \item{AccNo.
  CIDH920105}{Normalized average hydrophobicity scales (Cid
  et al., 1992)} \item{AccNo. BHAR880101}{Average
  flexibility indices (Bhaskaran-Ponnuswamy, 1988)}
  \item{AccNo. CHAM820101}{Polarizability parameter
  (Charton-Charton, 1982)} \item{AccNo. CHAM820102}{Free
  energy of solution in water, kcal/mole (Charton-Charton,
  1982)} \item{AccNo. CHOC760101}{Residue accessible
  surface area in tripeptide (Chothia, 1976)} \item{AccNo.
  BIGC670101}{Residue volume (Bigelow, 1967)} \item{AccNo.
  CHAM810101}{Steric parameter (Charton, 1981)}
  \item{AccNo. DAYM780201}{Relative mutability (Dayhoff et
  al., 1978b)}}}

  \item{nlag}{Maximum value of the lag parameter. Default
  is \code{30}.}
}
\value{
  A length \code{nlag} named vector
}
\description{
  Geary Autocorrelation Descriptor
}
\details{
  This function calculates the Geary autocorrelation
  descriptor (Dim: \code{length(props) * nlag}).
}
\examples{
A06852 = readFASTA(system.file('AAseq/A06852.fasta', package = 'rdpi'))
x = 'MPRLFSYLLGVWLLLSQLPREIPGQSTNDFIKACGRELVRLWVEICGSVSWGRTALSLEEPQLETGPPAETMPSSITKDAEILKMMLEFVPNLPQELKATLSERQPSLRELQQSASKDSNLNFEEFKKIILNRQNEAEDKSLLELKNLGLDKHSRKKRLFRMTLSEKCCQVGCIRKDIARLC'
extractGeary(x)
}
\author{
  Xiao Nan <\url{http://www.road2stat.com}>
}
\references{
  References needed here!
}
\seealso{
  See \code{\link{extractMoreauBroto}} and
  \code{\link{extractMoran}} for Moreau-Broto
  autocorrelation descriptors and Moran autocorrelation
  descriptors.
}
\keyword{autocorrelation}
\keyword{extract}
\keyword{extractGeary}
\keyword{Geary}

