#' Read Protein Sequences in PDB Format
#'
#' This function reads protein sequences in PDB (Protein Data Bank) format,
#' and return the amino acid sequences represented by single-letter code.
#'
#' @param file Path to the file containing the protein sequences in PDB format.
#' If it does not contain an absolute or
#' relative path, the file name is relative to the current
#' working directory, \code{\link{getwd}}.
#' The default here is to read the \code{4HHB.PDB} file which
#' is present in the \code{protseq} directory of the protr package.
#'
#' @return Character vector of the protein sequence.
#'
#' @keywords read PDB
#'
#' @aliases readPDB
#'
#' @author Nan Xiao <\url{https://nanx.me}>
#'
#' @seealso See \code{\link{readFASTA}} for reading
#' protein sequences in FASTA format.
#'
#' @export readPDB
#'
#' @references
#' Protein Data Bank Contents Guide: Atomic Coordinate Entry Format Description,
#' Version 3.30. Accessed 2013-06-26.
#' \url{ftp://ftp.wwpdb.org/pub/pdb/doc/format_descriptions/Format_v33_Letter.pdf}
#'
#' @examples
#' Seq4HHB = readPDB(system.file("protseq/4HHB.pdb", package = "protr"))

readPDB = function(
  file = system.file('protseq/4HHB.pdb', package = 'protr')) {

  pdb = scan(file, sep = '\n', what = 'complex', quiet = TRUE)

  category = gsub(' ', '', unlist(lapply(pdb, substring, 1, 6)))

  atomt  = unlist(lapply(pdb[category == 'ATOM'], substring, 14, 16))
  aminoa = as.character(unlist(lapply(pdb[category == 'ATOM'], substring, 18, 20)))
  chain  = unlist(lapply(pdb[category == 'ATOM'], substring, 22, 22))
  namino = as.numeric(unlist(lapply(pdb[category == 'ATOM'], substring, 23, 26)))
  atom   = data.frame(atom = atomt, aa = aminoa, chain = chain, naa = namino)

  uch = unique(chain)
  n   = length(uch)
  aaseq = vector('list', n)  # each element in the list is a chain

  for (i in 1:n) {
    subatom = atom[atom$chain == uch[i], ]
    uaa = unique(subatom$naa)
    aaseq[[i]] = as.factor(as.character(
      subatom$aa[unlist(lapply(
        uaa, FUN = function(X) { grep(X, subatom$naa)[1] } ))]))
    names(aaseq)[i] = paste('chain_', uch[i], sep = '')
  }

  dict = c(
    GLY = 'G', PRO = 'P', ALA = 'A', VAL = 'V', LEU = 'L',
    ILE = 'I', MET = 'M', CYS = 'C', PHE = 'F', TYR = 'Y',
    TRP = 'W', HIS = 'H', LYS = 'K', ARG = 'R', GLN = 'Q',
    ASN = 'N', GLU = 'E', ASP = 'D', SER = 'S', THR = 'T')

  for (i in 1:n) {
    pos = match(levels(aaseq[[i]]), names(dict))
    levels(aaseq[[i]]) = dict[pos]
    aaseq[[i]] = paste(as.character(aaseq[[i]]), collapse = '')
  }

  aaseq

}
