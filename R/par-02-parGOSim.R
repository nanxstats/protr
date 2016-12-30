.goPairSim = function (twoid, golist = golist, ont = ont, organism = organism, measure = measure, combine = combine) {

  id1 = twoid[1]
  id2 = twoid[2]

  if ( all(golist[[id1]] == '') | all(golist[[id2]] == '') ) {

    sim = 0L

  } else {

    id1good = 1:length(golist[[id1]])
    id2good = 1:length(golist[[id2]])

    gid1 = as.character(golist[[id1]][id1good])
    gid2 = as.character(golist[[id2]][id2good])

    res = try(suppressWarnings(GOSemSim::mgoSim(gid1, gid2,
                                                ont = ont,
                                                organism = organism,
                                                measure = measure,
                                                combine = combine)),
              silent = TRUE)

    if ( is.numeric(res) ) {
      sim = res
    } else {
      sim = 0L
    }

  }

  return(sim)

}

#' Protein Sequence Similarity Calculation based on Gene Ontology (GO) Similarity
#'
#' Protein Sequence Similarity Calculation based on Gene Ontology (GO) Similarity
#'
#' This function calculates
#' protein sequence similarity based on Gene Ontology (GO) similarity.
#'
#' @param golist A character vector, each component contains
#' a character vector of GO terms or one Entrez Gene ID.
#' @param type Input type of \code{golist}, \code{'go'} for GO Terms, \code{'gene'} for gene ID.
#' @param ont Default is \code{'MF'}, could be one of \code{'MF'}, \code{'BP'}, or \code{'CC'} subontologies.
#' @param organism Default is \code{'human'}, could be one of \code{'anopheles'}, \code{'arabidopsis'}, \code{'bovine'}, \code{'canine'},
#' \code{'chicken'}, \code{'chimp'}, \code{'coelicolor'}, \code{'ecolik12'},
#' \code{'ecsakai'}, \code{'fly'}, \code{'human'}, \code{'malaria'},
#' \code{'mouse'}, \code{'pig'}, \code{'rat'}, \code{'rhesus'},
#' \code{'worm'}, \code{'xenopus'}, \code{'yeast'} or \code{'zebrafish'}.
#' @param measure Default is \code{'Resnik'}, could be one of \code{'Resnik'}, \code{'Lin'}, \code{'Rel'}, \code{'Jiang'} or \code{'Wang'}.
#' @param combine Default is \code{'BMA'}, could be one of \code{'max'}, \code{'average'}, \code{'rcmax'} or \code{'BMA'}
#' for combining semantic similarity scores of multiple GO terms associated with protein.
#' @return A \code{n} x \code{n} similarity matrix.
#'
#' @keywords GO Ontology similarity
#'
#' @aliases parGOSim
#'
#' @author Nan Xiao <\url{http://nanx.me}>
#'
#' @seealso See \code{\link{twoGOSim}} for calculating the
#' GO semantic similarity between two groups of GO terms or two Entrez gene IDs.
#' See \code{\link{parSeqSim}} for paralleled protein similarity
#' calculation based on Smith-Waterman local alignment.
#'
#' @export parGOSim
#'
#' @examples
#' \dontrun{
#' # Be careful when testing this since it involves GO similarity computation
#' # and might produce unpredictable results in some environments
#'
#' require(GOSemSim)
#' require(org.Hs.eg.db)
#'
#' # by GO Terms
#' go1 = c('GO:0005215', 'GO:0005488', 'GO:0005515', 'GO:0005625', 'GO:0005802', 'GO:0005905')  # AP4B1
#' go2 = c('GO:0005515', 'GO:0005634', 'GO:0005681', 'GO:0008380', 'GO:0031202')  # BCAS2
#' go3 = c('GO:0003735', 'GO:0005622', 'GO:0005840', 'GO:0006412')  # PDE4DIP
#' glist = list(go1, go2, go3)
#' gsimmat1 = parGOSim(glist, type = 'go', ont = 'CC')
#' print(gsimmat1)
#'
#' # by Entrez gene id
#' genelist = list(c('150', '151', '152', '1814', '1815', '1816'))
#' gsimmat2 = parGOSim(genelist, type = 'gene')
#' print(gsimmat2)}

parGOSim = function (golist, type = c('go', 'gene'),
                     ont = 'MF', organism = 'human',
                     measure = 'Resnik', combine = 'BMA') {

  if ( type == 'gene' ) {
    gosimmat = GOSemSim::mgeneSim(unlist(golist), ont = ont, organism = organism, measure = measure, combine = combine, verbose = FALSE)
  }

  if ( type == 'go' ) {

    # generate lower matrix index
    idx = combn(1:length(golist), 2)

    # input is all pair combination
    gosimlist = vector('list', ncol(idx))

    for ( i in 1:ncol(idx) ) {
      gosimlist[[i]] = .goPairSim(rev(idx[, i]), golist = golist, ont = ont, organism = organism, measure = measure, combine = combine)
    }

    # convert list to matrix
    gosimmat = matrix(0, length(golist), length(golist))
    for (i in 1:length(gosimlist)) gosimmat[idx[2, i], idx[1, i]] = gosimlist[[i]]
    gosimmat[upper.tri(gosimmat)] = t(gosimmat)[upper.tri(t(gosimmat))]
    diag(gosimmat) = 1

  }

  return(gosimmat)

}

#' Protein Similarity Calculation based on Gene Ontology (GO) Similarity
#'
#' Protein Similarity Calculation based on Gene Ontology (GO) Similarity
#'
#' This function calculates the Gene Ontology (GO) similarity
#' between two groups of GO terms or two Entrez gene IDs.
#'
#' @param id1 A character vector. length > 1: each element is a GO term;
#' length = 1: the Entrez Gene ID.
#' @param id2 A character vector. length > 1: each element is a GO term;
#' length = 1: the Entrez Gene ID.
#' @param type Input type of id1 and id2, \code{'go'} for GO Terms, \code{'gene'} for gene ID.
#' @param ont Default is \code{'MF'}, could be one of \code{'MF'}, \code{'BP'}, or \code{'CC'} subontologies.
#' @param organism Default is \code{'human'}, could be one of \code{'anopheles'}, \code{'arabidopsis'}, \code{'bovine'}, \code{'canine'},
#' \code{'chicken'}, \code{'chimp'}, \code{'coelicolor'}, \code{'ecolik12'},
#' \code{'ecsakai'}, \code{'fly'}, \code{'human'}, \code{'malaria'},
#' \code{'mouse'}, \code{'pig'}, \code{'rat'}, \code{'rhesus'},
#' \code{'worm'}, \code{'xenopus'}, \code{'yeast'} or \code{'zebrafish'}.
#' @param measure Default is \code{'Resnik'}, could be one of \code{'Resnik'}, \code{'Lin'}, \code{'Rel'}, \code{'Jiang'} or \code{'Wang'}.
#' @param combine Default is \code{'BMA'}, could be one of \code{'max'}, \code{'average'}, \code{'rcmax'} or \code{'BMA'}
#' for combining semantic similarity scores of multiple GO terms associated with protein.
#' @return A n x n matrix.
#'
#' @keywords GO Ontology similarity
#'
#' @aliases twoGOSim
#'
#' @author Nan Xiao <\url{http://nanx.me}>
#'
#' @seealso See \code{\link{parGOSim}} for
#' protein similarity calculation based on
#' Gene Ontology (GO) semantic similarity.
#' See \code{\link{parSeqSim}} for paralleled protein similarity
#' calculation based on Smith-Waterman local alignment.
#'
#' @export twoGOSim
#'
#' @examples
#' \dontrun{
#' # Be careful when testing this since it involves GO similarity computation
#' # and might produce unpredictable results in some environments
#'
#' require(GOSemSim)
#' require(org.Hs.eg.db)
#'
#' # by GO terms
#' go1 = c("GO:0004022", "GO:0004024", "GO:0004023")
#' go2 = c("GO:0009055", "GO:0020037")
#' gsim1 = twoGOSim(go1, go2, type = 'go', ont = 'MF', measure = 'Wang')
#' print(gsim1)
#'
#' # by Entrez gene id
#' gene1 = '241'
#' gene2 = '251'
#' gsim2 = twoGOSim(gene1, gene2, type = 'gene', ont = 'BP', measure = 'Lin')
#' print(gsim2)}

twoGOSim = function (id1, id2, type = c('go', 'gene'),
                     ont = 'MF', organism = 'human',
                     measure = 'Resnik', combine = 'BMA') {

  if ( type == 'go' ) {
    sim = GOSemSim::mgoSim(id1, id2,
                           ont = ont, organism = organism,
                           measure = measure, combine = combine)
  }

  if ( type == 'gene' ) {
    sim = GOSemSim::geneSim(id1, id2,
                            ont = ont, organism = organism,
                            measure = measure, combine = combine)$geneSim
  }

  return(sim)

}
