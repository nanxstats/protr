\name{parGOSim}
\alias{parGOSim}
\title{Parallellized Protein Sequence Similarity Calculation based on Gene Ontology (GO) Similarity}
\usage{
  parGOSim(golist, type = c("go", "gene"), cores = 2,
    ont = "MF", organism = "human", measure = "Resnik",
    combine = "BMA")
}
\arguments{
  \item{golist}{A character vector, each component contains
  a character vector of GO terms or one Entrez Gene ID.}

  \item{type}{Input type of \code{golist}, \code{'go'} for
  GO Terms, \code{'gene'} for gene ID.}

  \item{ont}{Default is \code{'MF'}, could be one of
  \code{'MF'}, \code{'BP'}, or \code{'CC'} subontologies.}

  \item{organism}{Default is \code{'human'}, could be one
  of \code{'anopheles'}, \code{'arabidopsis'},
  \code{'bovine'}, \code{'canine'}, \code{'chicken'},
  \code{'chimp'}, \code{'coelicolor'}, \code{'ecolik12'},
  \code{'ecsakai'}, \code{'fly'}, \code{'human'},
  \code{'malaria'}, \code{'mouse'}, \code{'pig'},
  \code{'rat'}, \code{'rhesus'}, \code{'worm'},
  \code{'xenopus'}, \code{'yeast'} or \code{'zebrafish'}.}

  \item{measure}{Default is \code{'Resnik'}, could be one
  of \code{'Resnik'}, \code{'Lin'}, \code{'Rel'},
  \code{'Jiang'} or \code{'Wang'}.}

  \item{combine}{Default is \code{'BMA'}, could be one of
  \code{'max'}, \code{'average'}, \code{'rcmax'} or
  \code{'BMA'} for combining semantic similarity scores of
  multiple GO terms associated with protein.}
}
\value{
  A \code{n} x \code{n} similarity matrix.
}
\description{
  Parallellized Protein Sequence Similarity Calculation
  based on Gene Ontology (GO) Similarity
}
\details{
  This function implemented the parallellized version for
  calculating protein sequence similarity based on Gene
  Ontology (GO) similarity.
}
\note{
  Note
}
\examples{
\dontrun{
# by GO Terms
go1 = c("GO:0004022", "GO:0004024", "GO:0004023")
go2 = c("GO:0009055", "GO:0020037")
go3 = c("GO:0020039", "GO:0009056")
glist = list(go1, go2, go3)
gsimmat1 = parGOSim(glist, type = 'go', cores = 2)
print(gsimmat1)
# by Entrez gene id
genelist = list("835", "5261", "241")
gsimmat2 = parGOSim(genelist, type = 'gene', cores = 2)
print(gsimmat2) }
}
\author{
  Nan Xiao <\url{http://www.road2stat.com}>
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

