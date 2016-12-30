#' Read Protein Sequences in FASTA Format
#'
#' Read Protein Sequences in FASTA Format
#'
#' This function reads protein sequences in FASTA format.
#'
#' @param file The name of the file which the sequences in fasta format are
#'             to be read from. If it does not contain an absolute or
#'             relative path, the file name is relative to the current
#'             working directory, \code{\link{getwd}}.
#'             The default here is to read the \code{P00750.fasta} file which
#'             is present in the \code{protseq} directory of the protr package.
#'
#' @param legacy.mode If set to \code{TRUE}, lines starting with a semicolon ';'
#'                    are ignored. Default value is \code{TRUE}.
#' @param seqonly If set to \code{TRUE}, only sequences as returned without
#'                attempt to modify them or to get their names and
#'                annotations (execution time is divided approximately
#'                by a factor 3). Default value is \code{FALSE}.
#'
#' @return The result character vector
#'
#' The three returned argument are just different forms of the same output.
#' If one is interested in a Mahalanobis metric over the original data space,
#' the first argument is all she/he needs. If a transformation into another
#' space (where one can use the Euclidean metric) is preferred, the second
#' returned argument is sufficient. Using A and B is equivalent in the
#' following sense.
#'
#' @keywords read FASTA
#'
#' @aliases readFASTA
#'
#' @note Note that any different sets of instances (chunklets),
#'       e.g. {1, 3, 7} and {4, 6}, might belong to the
#'       same class and might belong to different classes.
#'
#' @author Nan Xiao <\url{http://nanx.me}>
#'
#' @seealso See \code{\link{getUniProt}} for retrieving protein sequences
#'          from uniprot.org
#'
#' @export readFASTA
#'
#' @references
#' Pearson, W.R. and Lipman, D.J. (1988)
#' Improved tools for biological sequence comparison.
#' \emph{Proceedings of the National Academy of Sciences
#' of the United States of America}, \bold{85}: 2444-2448
#'
#' @examples
#' P00750 = readFASTA(system.file('protseq/P00750.fasta', package = 'protr'))

readFASTA = function (file = system.file('protseq/P00750.fasta', package = 'protr'),
                      legacy.mode = TRUE, seqonly = FALSE) {

  # Read the fasta file as a vector of strings

  lines = readLines(file)

  # Remove comment lines starting with a semicolon ';'

  if (legacy.mode) {
    comments = grep("^;", lines)
    if (length(comments) > 0) {
      lines = lines[-comments]
    }
  }

  # Get the line numbers where sequences names are

  ind = which(substr(lines, 1L, 1L) == ">")

  # Compute the total number of sequences

  nseq = length(ind)

  if (nseq == 0) stop("no line starting with a > character found")

  # Localize sequence data

  start = ind + 1
  end = ind - 1
  end = c(end[-1], length(lines))

  # Read in sequences

  sequences = lapply(seq_len(nseq),
                     function(i) paste(lines[start[i]:end[i]], collapse = ""))

  if (seqonly) return(sequences)

  # Read in sequence names

  nomseq = lapply(seq_len(nseq), function (i) {
    firstword = strsplit(lines[ind[i]], " ")[[1]][1]
    substr(firstword, 2, nchar(firstword))
  })

  # Give the sequences names to the list elements

  names(sequences) = nomseq
  return(sequences)

}
