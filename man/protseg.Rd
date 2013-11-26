\name{protseg}
\alias{protseg}
\title{Protein Sequence Segmentation}
\usage{
  protseg(x,
    aa = c("A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"),
    k = 7)
}
\arguments{
  \item{x}{A character vector, as the input protein
  sequence.}

  \item{aa}{A character, the amino acid type. one of
  \code{'A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I',
  'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'}}

  \item{k}{A positive integer, specifys the window size
  (half of the window), default is 7.}
}
\value{
  A named list, each component contains one of the
  segmentations (a character string), names of the list
  components are the positions of the specified amino acid
  in the sequence.
}
\description{
  Protein Sequence Segmentation
}
\details{
  This function extracts the segmentations from the protein
  sequence.
}
\examples{
x = readFASTA(system.file('protseq/P00750.fasta', package = 'protr'))[[1]]
protseg(x, aa = 'R', k = 5)
}
\author{
  Nan Xiao <\url{http://www.road2stat.com}>
}
\keyword{acid}
\keyword{amino}
\keyword{protein}
\keyword{seg}
\keyword{segment}
\keyword{segmentation}
\keyword{sequence}

