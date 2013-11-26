\name{getUniProt}
\alias{getUniProt}
\title{Get Protein Sequences from UniProt by Protein ID}
\usage{
  getUniProt(id)
}
\arguments{
  \item{id}{A character vector, as the protein ID(s).}
}
\value{
  A list, each component contains one of the protein
  sequences.
}
\description{
  Get Protein Sequences from UniProt by Protein ID
}
\details{
  This function get protein sequences from uniprot.org by
  protein ID(s).
}
\examples{
ids = c('P00750', 'P00751', 'P00752')
\dontrun{getUniProt(ids)}
}
\author{
  Nan Xiao <\url{http://www.road2stat.com}>
}
\references{
  UniProt. \url{http://www.uniprot.org/}
}
\seealso{
  See \code{\link{readFASTA}} for reading FASTA format
  files.
}
\keyword{UniProt}
\keyword{getUniProt}

