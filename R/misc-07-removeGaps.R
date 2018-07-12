#' Remove or replace gaps from protein sequences.
#'
#' Remove/replace gaps or any irregular characters from protein sequences,
#' to make them suitable for feature extraction or sequence alignment
#' based similarity computation.
#'
#' @param x character vector, containing the input protein sequence(s).
#' @param pattern character string contains the gap (or other irregular)
#' character to be removed or replaced. Default is \code{"-"}.
#' For advanced usage, see \code{\link{gsub}}.
#' @param replacement a replacement for matched characters.
#' Default is \code{""} (remove the matched character).
#' @param ... addtional parameters for \code{\link{gsub}}.
#'
#' @return a vector of protein sequence(s) with gaps or
#' irregular characters removed/replaced.
#'
#' @author Nan Xiao <\url{https://nanx.me}>
#'
#' @export removeGaps
#'
#' @examples
#' # amino acid sequences that contain gaps ("-")
#' aaseq <- list(
#'   "MHGDTPTLHEYMLDLQPETTDLYCYEQLSDSSE-EEDEIDGPAGQAEPDRAHYNIVTFCCKCDSTLRLCVQS",
#'   "MHGDTPTLHEYMLDLQPETTDLYCYEQLNDSSE-EEDEIDGPAGQAEPDRAHYNIVTFCCKCDSTLRLCVQS"
#' )
#'
#' \dontrun{
#' #' # gaps create issues for alignment
#' parSeqSim(aaseq, cores = 1)
#'
#' # remove the gaps
#' nogapseq <- removeGaps(aaseq)
#' parSeqSim(nogapseq, cores = 1)}

removeGaps <- function(x, pattern = "-", replacement = "", ...)
  gsub(pattern, replacement, x, ...)
