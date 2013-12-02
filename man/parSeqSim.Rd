\name{parSeqSim}
\alias{parSeqSim}
\title{Parallellized Protein Sequence Similarity Calculation based on Sequence Alignment}
\usage{
  parSeqSim(protlist, cores = 2, type = "local",
    submat = "BLOSUM62")
}
\arguments{
  \item{protlist}{A length \code{n} list containing
  \code{n} protein sequences, each component of the list is
  a character string, storing one protein sequence. Unknown
  sequences should be represented as \code{''}.}

  \item{type}{Type of alignment, default is \code{'local'},
  could be \code{'global'} or \code{'local'}, where
  \code{'global'} represents Needleman-Wunsch global
  alignment; \code{'local'} represents Smith-Waterman local
  alignment.}

  \item{submat}{Substitution matrix, default is
  \code{'BLOSUM62'}, could be one of \code{'BLOSUM45'},
  \code{'BLOSUM50'}, \code{'BLOSUM62'}, \code{'BLOSUM80'},
  \code{'BLOSUM100'}, \code{'PAM30'}, \code{'PAM40'},
  \code{'PAM70'}, \code{'PAM120'}, \code{'PAM250'}.}
}
\value{
  A \code{n} x \code{n} similarity matrix.
}
\description{
  Parallellized Protein Sequence Similarity Calculation
  based on Sequence Alignment
}
\details{
  This function implemented the parallellized version for
  calculating protein sequence similarity based on sequence
  alignment.
}
\note{
  Note
}
\examples{
\dontrun{
s1 = readFASTA(system.file('protseq/P00750.fasta', package = 'protr'))[[1]]
s2 = readFASTA(system.file('protseq/P08218.fasta', package = 'protr'))[[1]]
s3 = readFASTA(system.file('protseq/P10323.fasta', package = 'protr'))[[1]]
s4 = readFASTA(system.file('protseq/P20160.fasta', package = 'protr'))[[1]]
s5 = readFASTA(system.file('protseq/Q9NZP8.fasta', package = 'protr'))[[1]]
plist = list(s1, s2, s3, s4, s5)
psimmat = parSeqSim(plist, cores = 2, type = 'local', submat = 'BLOSUM62')
print(psimmat) }
}
\author{
  Nan Xiao <\url{http://www.road2stat.com}>
}
\seealso{
  See \code{\link{parGOSim}} for paralleled protein
  similarity calculation based on Gene Ontology (GO)
  similarity.
}
\keyword{Needleman-Wunsch}
\keyword{Smith-Waterman}
\keyword{alignment}
\keyword{global}
\keyword{local}
\keyword{parSeqSim}
\keyword{parallel}
\keyword{sequence}
\keyword{similarity}

