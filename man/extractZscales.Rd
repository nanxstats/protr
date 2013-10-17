\name{extractZscales}
\alias{extractZscales}
\title{Z-scales Descriptor}
\usage{
  extractZscales(x, pc = 3, lag = 5, scale = TRUE,
    customprops = NULL)
}
\arguments{
  \item{x}{A character vector, as the input protein
  sequence.}
}
\value{
  A length \code{20 + n * lambda} named vector, \code{n} is
  the number of properties selected.
}
\description{
  Z-scales Descriptor
}
\details{
  This function calculates the Z-scales descriptor (Dim:
  \code{20 + (n * lambda)}, \code{n} is the number of
  properties selected, default is 80).
}
\note{
  Note
}
\examples{
x = readFASTA(system.file('protseq/P00750.fasta', package = 'protr'))[[1]]
extractZscales(x)
}
\author{
  Xiao Nan <\url{http://www.road2stat.com}>
}
\references{
  Sandberg, M., Eriksson, L., Jonsson, J.,
  Sj\"{o}str\"{o}m, M., & Wold, S. (1998). New chemical
  descriptors relevant for the design of biologically
  active peptides. A multivariate characterization of 87
  amino acids. \emph{Journal of medicinal chemistry},
  41(14), 2481--2491.

  Hellberg, S., Sj\"{o}str\"{o}m, M., Skagerberg, B., &
  Wold, S. (1987). Peptide quantitative structure-activity
  relationships, a multivariate approach. \emph{Journal of
  medicinal chemistry}, 30(7), 1126--1135.
}
\seealso{
  See \code{\link{extractTscales}},
  \code{\link{extractSTscales}} and
  \code{\link{extractVHSE}} for T-scales, ST-scales and
  VHSE descriptors. See \code{\link{extractScales}} for
  generalized scales-based descriptors.
}
\keyword{extract}
\keyword{extractZscales}
\keyword{PCM}
\keyword{Z-scales}

