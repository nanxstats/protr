#' Get Protein Sequences from UniProt by Protein ID
#'
#' Get Protein Sequences from UniProt by Protein ID
#' 
#' This function get protein sequences from uniprot.org by protein ID(s).
#' 
#' @param id A character vector, as the protein ID(s). 
#'
#' @return A list, each component contains one of the protein sequences. 
#' 
#' @keywords UniProt getUniProt
#'
#' @aliases getUniProt
#' 
#' @author Xiao Nan <\url{http://www.road2stat.com}>
#' 
#' @seealso See \code{\link{readFASTA}} for reading FASTA format files.
#' 
#' @export getUniProt
#' 
#' @references
#' UniProt. \url{http://www.uniprot.org/}
#' 
#' @examples
#' ids = c('P00750', 'P00751', 'P00752')
#' \dontrun{getUniProt(ids)}
#' 

getUniProt = function (id) {
  
  id = as.character(id)
  
  n = length(id)
  
  proteins = vector('list', n)
  
  # API format:
  # http://www.uniprot.org/uniprot/P00750.fasta
  
  for (i in 1:n) {
    proteins[[i]] = readFASTA(paste('http://www.uniprot.org/uniprot/', 
                                    id[i], '.fasta', sep = ''))[[1]]
  }
    
  return(proteins)
  
}

