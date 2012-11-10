\name{extractAPAAC}
\alias{extractAPAAC}
\title{Amphiphilic Pseudo Amino Acid Composition Descriptor}
\usage{
  extractAPAAC(x,
    props = c("Hydrophobicity", "Hydrophilicity"),
    lambda = 30, w = 0.05)
}
\arguments{
  \item{x}{A character vector, as the input protein
  sequence.}

  \item{props}{A character vector, specifying the
  properties used.  2 properties are used by default, as
  listed below: \describe{
  \item{\code{'Hydrophobicity'}}{Hydrophobicity value of
  the 20 amino acids}
  \item{\code{'Hydrophilicity'}}{Hydrophilicity value of
  the 20 amino acids}}}

  \item{lambda}{The lambda parameter for the APAAC
  descriptors, default is 30.}

  \item{w}{The weighting factor, default is 0.05.}
}
\value{
  A length \code{20 + n * lambda} named vector, \code{n} is
  the number of properties selected.
}
\description{
  Amphiphilic Pseudo Amino Acid Composition Descriptor
}
\details{
  This function calculates the Amphiphilic Pseudo Amino
  Acid Composition (APAAC) descriptor (Dim: \code{20 + (n *
  lambda)}, \code{n} is the number of properties selected,
  default is 80).
}
\note{
  Note the default \code{20 * 2} \code{prop} values have
  been already independently given in the function. Users
  could also specify other (up to 544) properties with the
  Accession Number in the \code{\link{AAindex}} data, with
  or without the default three properties, which means
  users should explicitly specify the properties to use.
}
\examples{
x = readFASTA(system.file('AAseq/P00750.fasta', package = 'rdpi'))[[1]]
extractAPAAC(x)
}
\author{
  Xiao Nan <\url{http://www.road2stat.com}>
}
\references{
  Kuo-Chen Chou. Prediction of Protein Cellular Attributes
  Using Pseudo-Amino Acid Composition. \emph{PROTEINS:
  Structure, Function, and Genetics}, 2001, 43: 246-255.

  Type 2 pseudo amino acid composition.
  \url{http://www.csbio.sjtu.edu.cn/bioinf/PseAAC/type2.htm}

  Kuo-Chen Chou. Using Amphiphilic Pseudo Amino Acid
  Composition to Predict Enzyme Subfamily Classes.
  \emph{Bioinformatics}, 2005, 21, 10-19.

  JACS, 1962, 84: 4240-4246. (C. Tanford). (The
  hydrophobicity data)

  PNAS, 1981, 78:3824-3828 (T.P.Hopp & K.R.Woods). (The
  hydrophilicity data)
}
\seealso{
  See \code{\link{extractPAAC}} for pseudo amino acid
  composition descriptor.
}
\keyword{Acid}
\keyword{Amino}
\keyword{Amphiphilic}
\keyword{APAAC}
\keyword{Composition}
\keyword{extract}
\keyword{extractAPAAC}
\keyword{Pseudo}

