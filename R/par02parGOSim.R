.goPairSim = function (twoid, ont = ont, organism = organism, measure = measure, combine = combine) {
  
  id1 = twoid[1]
  id2 = twoid[2]
  
  if ( all(golist[[id1]] == '') | all(golist[[id2]] == '') ) {
    
    sim = 0L
    
  } else {
    
    id1good = 1:length(golist[[id1]])
    id2good = 1:length(golist[[id2]])
    
    gid1 = as.character(golist[[id1]][id1good])
    gid2 = as.character(golist[[id2]][id2good])
        
    res = try(suppressWarnings(mgoSim(gid1, gid2, ont = ont, organism = organism, measure = measure, combine = combine)), silent = TRUE)
    
    if ( is.list(res) ) {
      if ( res$geneSim == Inf ) {
        sim = 1L
        } else {
          sim = res$geneSim
        }
      } else {
        sim = 0L
      }
    
  }
  
  return(sim)
  
}

#' Parallellized Protein Sequence Similarity Calculation based on Gene Ontology (GO) Similarity
#'
#' Parallellized Protein Sequence Similarity Calculation based on Gene Ontology (GO) Similarity
#'
#' This function implemented the parallellized version for calculating 
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
#' @return A nxn matrix.
#' 
#' @note Note
#' 
#' @keywords GO Gene Ontology parallel similarity parGOSim
#'
#' @aliases parGOSim
#' 
#' @author Nan Xiao <\url{http://www.road2stat.com}>
#' 
#' @seealso See \code{\link{parSeqSim}} for paralleled protein similarity
#' calculation based on Smith-Waterman local alignment.
#' 
#' @export parGOSim
#' 
#' @examples
#' go1 = c("GO:0004022", "GO:0004024", "GO:0004023")
#' go2 = c("GO:0009055", "GO:0020037")
#' go3 = c("GO:0020039", "GO:0009056")
#' golist = list(go1, go2, go3)
#' genelist = list("835", "5261", "241")
#' parGOSim(golist, type = 'go', cores = 4, ont = 'MF', organism = 'human', measure = 'Resnik', combine = 'BMA')
#' parGOSim(genelist, type = 'gene', ont = 'MF', organism = 'human', measure = 'Resnik', combine = 'BMA')
#' 

parGOSim = function (golist, type = c('go', 'gene'), 
                     cores = 2, 
                     ont = 'MF', organism = 'human', 
                     measure = 'Resnik', combine = 'BMA') {
  
  require(GOSemSim)
  require(foreach)
  if (.Platform$OS.type == 'windows') {
    require(doParallel)
    registerDoParallel(cores)
  } else {
    require(doMC)
    registerDoMC(cores)
  }
  
  if ( type == 'gene' ) {
    gosimmat = mgeneSim(unlist(golist), ont = ont, organism = organism, measure = measure, combine = combine, verbose = FALSE)
  }
  
  if ( type == 'go' ) {
    
    # generate lower matrix index
    idx = combn(1:length(golist), 2)
    
    # then use foreach parallelization
    # input is all pair combination
    
    gosimlist = vector('list', ncol(idx))
    
    gosimlist <- foreach (i = 1:ncol(idx), .errorhandling = 'pass') %dopar% {
      xxx <- .goPairSim(rev(idx[, i]), ont = ont, organism = organism, measure = measure, combine = combine)
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
#' @return A nxn matrix.
#' 
#' @note Note
#' 
#' @keywords GO Gene Ontology similarity twoGOSim
#'
#' @aliases twoGOSim
#' 
#' @author Nan Xiao <\url{http://www.road2stat.com}>
#' 
#' @seealso See \code{\link{parSeqSim}} for paralleled protein similarity
#' calculation based on Smith-Waterman local alignment.
#' 
#' @export twoGOSim
#' 
#' @examples
#' go1 = c("GO:0004022", "GO:0004024", "GO:0004023")
#' go2 = c("GO:0009055", "GO:0020037")
#' twoGOSim(go1, go2, type = 'go', ont = 'MF', organism = 'human', measure = 'Resnik', combine = 'BMA')
#' gene1 = '241'
#' gene2 = '251'
#' twoGOSim(gene1, gene2, type = 'gene', ont = 'MF', organism = 'human', measure = 'Resnik', combine = 'BMA')
#' 

twoGOSim = function (id1, id2, type = c('go', 'gene'), 
                     ont = 'MF', organism = 'human', 
                     measure = 'Resnik', combine = 'BMA') {
  
  if ( type == 'go' ) {
    sim = mgoSim(id1, id2, ont = ont, organism = organism, measure = measure)
  }
  
  if ( type == 'gene' ) {
    sim = geneSim(id1, id2, ont = ont, organism = organism, measure = measure)$geneSim
  }
  
  return(sim)
  
}
