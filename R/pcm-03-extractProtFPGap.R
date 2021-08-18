#' Amino Acid Properties Based Scales Descriptors (Protein Fingerprint)
#' with Gap Support
#'
#' This function calculates amino acid properties based scales descriptors
#' (protein fingerprint) with gap support. Users can specify which AAindex
#' properties to select from the AAindex database by specify the numerical or
#' character index of the properties in the AAindex database.
#'
#' @param x A character vector, as the input protein sequence.
#' Use '\code{-}' to represent gaps in the sequence.
#' @param index Integer vector or character vector. Specify which AAindex
#' properties to select from the AAindex database by specify the
#' numerical or character index of the properties in the AAindex database.
#' Default is \code{NULL}, means selecting all the AA properties
#' in the AAindex database.
#' @param pc Integer. Use the first pc principal components as the scales.
#' Must be no greater than the number of AA properties provided.
#' @param lag The lag parameter. Must be less than the amino acids.
#' @param scale Logical. Should we auto-scale the property matrix
#' before PCA? Default is \code{TRUE}.
#' @param silent Logical. Whether we print the standard deviation,
#' proportion of variance and the cumulative proportion of
#' the selected principal components or not. Default is \code{TRUE}.
#'
#' @return A length \code{lag * p^2} named vector,
#' \code{p} is the number of scales (principal components) selected.
#'
#' @author Nan Xiao <\url{https://nanx.me}>
#'
#' @importFrom stats na.omit
#'
#' @export extractProtFPGap
#'
#' @examples
#' # amino acid sequence with gaps
#' x <- readFASTA(system.file("protseq/align.fasta", package = "protr"))$`IXI_235`
#' fp <- extractProtFPGap(x, index = c(160:165, 258:296), pc = 5, lag = 7, silent = FALSE)
extractProtFPGap <- function(
  x, index = NULL, pc, lag, scale = TRUE, silent = TRUE) {
  propmat <- get("AAindex")

  if (!is.null(index)) {
    propmat <- t(na.omit(as.matrix(propmat[index, 7:26])))
  }

  res <- extractScalesGap(
    x = x, propmat = propmat, pc = pc,
    lag = lag, scale = scale, silent = silent
  )

  res
}
