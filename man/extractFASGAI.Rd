\name{extractFASGAI}
\alias{extractFASGAI}
\title{Factor Analysis Scales of Generalized AA Information (FASGAI) Descriptor}
\usage{
  extractFASGAI(x, propmat, factors, scores, lag, silent)
}
\arguments{
  \item{x}{A character vector, as the input protein
  sequence.}

  \item{propmat}{A matrix containing the properties for the
  amino acids.  Each row represent one amino acid type,
  each column represents one property.  Note that the
  one-letter row names must be provided for we need them to
  seek the properties for each AA type.}

  \item{factors}{Integer. The number of factors to be
  fitted.  Must be no greater than the number of AA
  properties provided.}

  \item{scores}{Type of scores to produce. The default is
  \code{"regression"}, which gives Thompson's scores,
  \code{"Bartlett"} given Bartlett's weighted least-squares
  scores.}

  \item{lag}{The lag parameter. Must be less than the amino
  acids.}

  \item{silent}{Logical. Whether we print the SS loadings,
  proportion of variance and the cumulative proportion of
  the selected factors or not.  Default is \code{TRUE}.}
}
\value{
  A length \code{lag * p^2} named vector, \code{p} is the
  number of scales (factors) selected.
}
\description{
  Factor Analysis Scales of Generalized AA Information
  (FASGAI) Descriptor
}
\details{
  This function calculates the (generalized) Factor
  Analysis Scales of Generalized AA Information (FASGAI)
  descriptors. Users could provide customized amino acid
  property matrices.
}
\note{
  Note
}
\examples{
x = readFASTA(system.file('protseq/P00750.fasta', package = 'protr'))[[1]]
# extractFASGAI(x)
}
\author{
  Nan Xiao <\url{http://www.road2stat.com}>
}
\references{
  Liang, G., & Li, Z. (2007). Factor analysis scale of
  generalized amino acid information as the source of a new
  set of descriptors for elucidating the structure and
  activity relationships of cationic antimicrobial
  peptides. QSAR & Combinatorial Science, 26(6), 754--763.
}
\keyword{FASGAI}
\keyword{PCM}
\keyword{extract}
\keyword{extractFASGAI}

