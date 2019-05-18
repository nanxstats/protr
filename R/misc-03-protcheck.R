#' Protein sequence amino acid type sanity check
#'
#' This function checks if the protein sequence's amino acid types
#' are in the 20 default types.
#'
#' @param x A character vector, as the input protein sequence.
#'
#' @return Logical. \code{TRUE} if all of the amino acid types
#' of the sequence are within the 20 default types.
#'
#' @author Nan Xiao <\url{https://nanx.me}>
#'
#' @export protcheck
#'
#' @examples
#' x <- readFASTA(system.file("protseq/P00750.fasta", package = "protr"))[[1]]
#' protcheck(x) # TRUE
#' protcheck(paste(x, "Z", sep = "")) # FALSE
protcheck <- function(x) {
  AADict <- c(
    "A", "R", "N", "D", "C", "E", "Q", "G", "H", "I",
    "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"
  )

  all(strsplit(x, split = "")[[1]] %in% AADict)
}
