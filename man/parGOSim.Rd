\name{parGOSim}
\alias{parGOSim}
\title{Large-Scale Parallelled Protein Sequence Similarity Calculation based on Gene Ontology (GO) Similarity}
\usage{
  parGOSim(x, propmat, pc, lag, scale = TRUE,
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

  \item{pc}{Integer. Use the first pc principal components
  as the scales.  Must be no greater than the number of AA
  properties provided.}

  \item{lag}{The lag parameter. Must be less than the amino
  acids.}

  \item{scale}{Logical. Should we auto-scale the property
  matrix (\code{propmat}) before PCA? Default is
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
  Large-Scale Parallelled Protein Sequence Similarity
  Calculation based on Gene Ontology (GO) Similarity
}
\details{
  This function implemented the parallelled version for
  calculating protein sequence similarity based on Gene
  Ontology (GO) similarity.
}
\note{
  Note
}
\examples{
x = readFASTA(system.file('protseq/P00750.fasta', package = 'protr'))[[1]]
AAidxmat = t(na.omit(as.matrix(AAindex[, 7:26])))
scales = parGOSim(x, propmat = AAidxmat, pc = 5, lag = 7, silent = FALSE)
}
\author{
  Xiao Nan <\url{http://www.road2stat.com}>
}
\seealso{
  See \code{\link{parSeqSim}} for paralleled protein
  similarity calculation based on Smith-Waterman local
  alignment.
}
\keyword{GO}
\keyword{Gene}
\keyword{Ontology}
\keyword{parGOSim}
\keyword{parallel}
\keyword{similarity}

