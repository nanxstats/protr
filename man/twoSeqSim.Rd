\name{twoSeqSim}
\alias{twoSeqSim}
\title{Protein Sequence Alignment for Two Protein Sequences}
\usage{
  twoSeqSim(seq1, seq2, type = "local",
    submat = "BLOSUM62")
}
\arguments{
  \item{seq1}{A character string, containing one protein
  sequence.}

  \item{seq2}{A character string, containing another
  protein sequence.}

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
  An Biostrings object containing the scores and other
  alignment information.
}
\description{
  Protein Sequence Alignment for Two Protein Sequences
}
\details{
  This function implements the sequence alignment between
  two protein sequences.
}
\note{
  Note
}
\examples{
\dontrun{
require(Biostrings)
twoSeqSim(seq1, seq2, type = 'local', submat = 'BLOSUM62') }
}
\author{
  Nan Xiao <\url{http://www.road2stat.com}>
}
\seealso{
  See \code{\link{parSeqSim}} for paralleled pairwise
  protein similarity calculation based on sequence
  alignment.
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

