#' Scales-Based Descriptors with 19 classes of Molecular Descriptors from DRAGON
#'
#' Scales-Based Descriptors with 19 classes of Molecular Descriptors from DRAGON
#'
#' This function calculates the scales-based descriptors with DRAGON descriptors.
#' Users could provide customized amino acid property matrices.
#' 
#' @param x A character vector, as the input protein sequence.
#' @param propmat A matrix containing the properties for the amino acids. 
#'        Each row represent one amino acid type, each column represents one property.
#'        Note that the one-letter row names must be provided for we need them to seek 
#'        the properties for each AA type.
#' @param k Integer. The maximum dimension of the space which the data 
#'        are to be represented in.
#'        Must be no greater than the number of AA properties provided.
#' @param lag The lag parameter. Must be less than the amino acids.
#' @param scale Logical. Should we auto-scale the property matrix 
#'        (\code{propmat}) before doing MDS? Default is \code{TRUE}.
#' @param silent Logical. Whether we print the standard deviation, 
#'        proportion of variance and the cumulative proportion of 
#'        the selected principal components or not.
#'        Default is \code{TRUE}.
#'
#' @return A length \code{lag * p^2} named vector, 
#'         \code{p} is the number of scales selected.
#' 
#' @note Note
#' 
#' @keywords extract scales extractDescScales PCM
#'
#' @aliases extractDescScales
#' 
#' @author Xiao Nan <\url{http://www.road2stat.com}>
#' 
#' @seealso See \code{\link{extractScales}} for generalized 
#' AA-descriptor based scales descriptors.
#' 
#' @export extractDescScales
#' 
#' @examples
#' x = readFASTA(system.file('protseq/P00750.fasta', package = 'protr'))[[1]]
#' # AAidxmat = t(na.omit(as.matrix(AAindex[, 7:26])))
#' # scales = extractDescScales(x, propmat = AAidxmat, pc = 5, lag = 7, silent = FALSE)
#' 

extractDescScales = function (x, propmat, k, lag, scale = TRUE, silent = TRUE) {
  
  # 提供一个列表(17 18 个数据) 来自于DRAGON计算的描述符 用户可以选择
  AAidxmat = t(na.omit(as.matrix(AAindex[, 7:26])))
  
  d = dist(AAidxmat) # euclidean distances between the rows
  fit = cmdscale(d, eig = TRUE, k = 5) # k is the number of dim
  fit
  
}
