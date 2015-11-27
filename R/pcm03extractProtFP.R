#' Amino Acid Properties Based Scales Descriptors (Protein Fingerprint)
#'
#' Amino Acid Properties Based Scales Descriptors (Protein Fingerprint)
#'
#' This function calculates amino acid properties based scales descriptors
#' (protein fingerprint). Users could specify which AAindex properties to
#' select from the AAindex database by specify the numerical or character
#' index of the properties in the AAindex database.
#'
#' @param x A character vector, as the input protein sequence.
#' @param index Integer vector or character vector. Specify which AAindex
#'        properties to select from the AAindex database by specify the
#'        numerical or character index of the properties in the
#'        AAindex database.
#'        Default is \code{NULL}, means selecting all the AA properties
#'        in the AAindex database.
#' @param pc Integer. Use the first pc principal components as the scales.
#'        Must be no greater than the number of AA properties provided.
#' @param lag The lag parameter. Must be less than the amino acids.
#' @param scale Logical. Should we auto-scale the property matrix
#'        before PCA? Default is \code{TRUE}.
#' @param silent Logical. Whether we print the standard deviation,
#'        proportion of variance and the cumulative proportion of
#'        the selected principal components or not.
#'        Default is \code{TRUE}.
#'
#' @return A length \code{lag * p^2} named vector,
#'         \code{p} is the number of scales (principal components) selected.
#'
#' @keywords extract scales AAindex
#'
#' @aliases extractProtFP
#'
#' @author Nan Xiao <\url{http://nanx.me}>
#'
#' @export extractProtFP
#'
#' @examples
#' x = readFASTA(system.file('protseq/P00750.fasta', package = 'protr'))[[1]]
#' fp = extractProtFP(x, index = c(160:165, 258:296), pc = 5, lag = 7, silent = FALSE)

extractProtFP = function (x, index = NULL, pc, lag, scale = TRUE, silent = TRUE) {

  AAindex = get('AAindex')

  if (!is.null(index)) propmat = t(na.omit(as.matrix(AAindex[index, 7:26])))

  result = extractScales(x = x, propmat = propmat, pc = pc, lag = lag, scale = scale, silent = silent)

  return(result)

}
