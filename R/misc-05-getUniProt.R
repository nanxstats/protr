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
#' @keywords UniProt
#'
#' @aliases getUniProt
#'
#' @author Nan Xiao <\url{http://nanx.me}>
#'
#' @seealso See \code{\link{readFASTA}} for reading FASTA format files.
#'
#' @export getUniProt
#'
#' @examples
#' \donttest{
#' # Network latency may slow down this example
#' # Only test this when your connection is fast enough
#' ids = c('P00750', 'P00751', 'P00752')
#' getUniProt(ids)}

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
