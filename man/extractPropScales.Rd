\name{extractPropScales}
\alias{extractPropScales}
\title{Generalized AA-Properties Based Scales Descriptors}
\usage{
  extractPropScales(x, descmat, pc, lag, scale = TRUE,
    silent = TRUE)
}
\arguments{
  \item{x}{A character vector, as the input protein
  sequence.}

  \item{descmat}{A matrix containing the properties for the
  amino acids.  Each row represent one amino acid type,
  each column represents one property.  Note that the
  one-letter row names must be provided for we need them to
  seek the properties for each AA type.}

  \item{pc}{Integer. Use the first pc principal components
  as the scales.  Must be no greater than the number of AA
  properties provided.}

  \item{lag}{The lag parameter. Must be less than the amino
  acids.}

  \item{scale}{Logical. Should we auto-scale the property
  matrix (\code{descmat}) before PCA? Default is
  \code{TRUE}.}

  \item{silent}{Logical. Whether we print the standard
  deviation, proportion of variance and the cumulative
  proportion of the selected principal components or not.
  Default is \code{TRUE}.}
}
\value{
  A length \code{lag * p^2} named vector, \code{p} is the
  number of scales (principal components) selected.
}
\description{
  Generalized AA-Properties Based Scales Descriptors
}
\details{
  This function calculates the generalized amino acid
  properties based scales descriptors. Users could provide
  customized amino acid descriptor matrices.
}
\note{
  Note
}
\examples{
x = readFASTA(system.file('protseq/P00750.fasta', package = 'protr'))[[1]]
# data(AAindex)
# AAidxmat = t(na.omit(as.matrix(AAindex[, 7:26])))
# scales = extractPropScales(x, descmat = AAidxmat, pc = 5, lag = 7, silent = FALSE)
}
\author{
  Xiao Nan <\url{http://www.road2stat.com}>
}
\seealso{
  See \code{\link{extractScales}} for generalized
  scales-based descriptors.
}
\keyword{PCM}
\keyword{descriptor}
\keyword{extract}
\keyword{extractPropScales}
\keyword{scales}

