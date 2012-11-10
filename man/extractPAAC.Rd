\name{extractPAAC}
\alias{extractPAAC}
\title{Pseudo Amino Acid Composition Descriptor}
\usage{
  extractPAAC(x,
    props = c("Hydrophobicity", "Hydrophilicity", "SideChainMass"),
    lambda = 30, w = 0.05)
}
\arguments{
  \item{x}{A character vector, as the input protein
  sequence.}

  \item{props}{A character vector, specifying the
  properties used.  3 properties are used by default, as
  listed below: \describe{
  \item{\code{'Hydrophobicity'}}{Hydrophobicity value of
  the 20 amino acids}
  \item{\code{'Hydrophilicity'}}{Hydrophilicity value of
  the 20 amino acids}
  \item{\code{'SideChainMass'}}{Side-chain mass of the 20
  amino acids}} Note these \code{20 * 3} values have been
  already independently given in the function.  Users could
  also specify other (up to 544) properties with the
  Accession Number in the \code{\link{AAindex}} data, with
  or without the above three properties, which means users
  should explicitly specify the properties to use.}

  \item{lambda}{The lambda parameter for the PAAC
  descriptors, default is 30.}

  \item{w}{The weighting factor, default is 0.05.}
}
\value{
  A length \code{20 + lambda} named vector
}
\description{
  Pseudo Amino Acid Composition Descriptor
}
\details{
  This function calculates the Pseudo Amino Acid
  Composition (PAAC) descriptor (Dim: \code{20 + lambda},
  default is 50).
}
\examples{
x = readFASTA(system.file('AAseq/P00750.fasta', package = 'rdpi'))[[1]]
extractPAAC(x)
}
\author{
  Xiao Nan <\url{http://www.road2stat.com}>
}
\references{
  Kuo-Chen Chou. Prediction of Protein Cellular Attributes
  Using Pseudo-Amino Acid Composition. \emph{PROTEINS:
  Structure, Function, and Genetics}, 2001, 43: 246-255.

  Type 1 pseudo amino acid composition.
  \url{http://www.csbio.sjtu.edu.cn/bioinf/PseAAC/type1.htm}

  Kuo-Chen Chou. Using Amphiphilic Pseudo Amino Acid
  Composition to Predict Enzyme Subfamily Classes.
  \emph{Bioinformatics}, 2005, 21, 10-19.

  JACS, 1962, 84: 4240-4246. (C. Tanford). (The
  hydrophobicity data)

  PNAS, 1981, 78:3824-3828 (T.P.Hopp & K.R.Woods). (The
  hydrophilicity data)

  CRC Handbook of Chemistry and Physics, 66th ed., CRC
  Press, Boca Raton, Florida (1985). (The side-chain mass
  data)

  R.M.C. Dawson, D.C. Elliott, W.H. Elliott, K.M. Jones,
  Data for Biochemical Research 3rd ed., Clarendon Press
  Oxford (1986). (The side-chain mass data)
}
\seealso{
  See \code{\link{extractAPAAC}} for amphiphilic pseudo
  amino acid composition descriptor.
}
\keyword{Acid}
\keyword{Amino}
\keyword{Composition}
\keyword{extract}
\keyword{extractPAAC}
\keyword{PAAC}
\keyword{Pseudo}

