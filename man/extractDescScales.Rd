\name{extractDescScales}
\alias{extractDescScales}
\title{Scales-Based Descriptors with 20+ classes of Molecular Descriptors from DRAGON, DS and MOE}
\usage{
  extractDescScales(x, propmat, k, lag, scale = TRUE,
    silent = TRUE)
}
\arguments{
  \item{x}{A character vector, as the input protein
  sequence.}

  \item{propmat}{A matrix containing the properties for the
  amino acids.  Each row represent one amino acid type,
  each column represents one property.  Note that the
  one-letter row names must be provided for we need them to
  seek the properties for each AA type.}

  \item{k}{Integer. The maximum dimension of the space
  which the data are to be represented in.  Must be no
  greater than the number of AA properties provided.}

  \item{lag}{The lag parameter. Must be less than the amino
  acids.}

  \item{scale}{Logical. Should we auto-scale the property
  matrix (\code{propmat}) before doing MDS? Default is
  \code{TRUE}.}

  \item{silent}{Logical. Whether we print the standard
  deviation, proportion of variance and the cumulative
  proportion of the selected principal components or not.
  Default is \code{TRUE}.}
}
\value{
  A length \code{lag * p^2} named vector, \code{p} is the
  number of scales selected.
}
\description{
  Scales-Based Descriptors with 20+ classes of Molecular
  Descriptors from DRAGON, DS and MOE
}
\details{
  This function calculates the scales-based descriptors
  with DRAGON descriptors. Users could provide customized
  amino acid property matrices.
}
\note{
  Note
}
\examples{
x = readFASTA(system.file('protseq/P00750.fasta', package = 'protr'))[[1]]
# AAidxmat = t(na.omit(as.matrix(AAindex[, 7:26])))
# scales = extractDescScales(x, propmat = AAidxmat, pc = 5, lag = 7, silent = FALSE)
}
\author{
  Nan Xiao <\url{http://www.road2stat.com}>
}
\seealso{
  See \code{\link{extractScales}} for generalized
  AA-descriptor based scales descriptors.
}
\keyword{PCM}
\keyword{extract}
\keyword{extractDescScales}
\keyword{scales}

