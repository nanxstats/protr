#' Factor Analysis Scales of Generalized AA Information (FASGAI) Descriptor
#'
#' Factor Analysis Scales of Generalized AA Information (FASGAI) Descriptor
#' 
#' This function calculates the (generalized) Factor Analysis Scales 
#' of Generalized AA Information (FASGAI) descriptors.
#' Users could provide customized amino acid property matrices.
#' 
#' @param x A character vector, as the input protein sequence.
#' @param propmat A matrix containing the properties for the amino acids. 
#'        Each row represent one amino acid type, each column represents one property.
#'        Note that the one-letter row names must be provided for we need them to seek 
#'        the properties for each AA type.
#' @param factors Integer. The number of factors to be fitted.
#'        Must be no greater than the number of AA properties provided.
#' @param scores Type of scores to produce. The default is \code{"regression"}, 
#'        which gives Thompson's scores, \code{"Bartlett"} given Bartlett's weighted 
#'        least-squares scores.
#' @param lag The lag parameter. Must be less than the amino acids.
#' @param silent Logical. Whether we print the SS loadings, 
#'        proportion of variance and the cumulative proportion of 
#'        the selected factors or not.
#'        Default is \code{TRUE}.
#'        
#' @return A length \code{lag * p^2} named vector, 
#'         \code{p} is the number of scales (factors) selected.
#' 
#' @note Note
#' 
#' @keywords extract FASGAI extractFASGAI PCM
#'
#' @aliases extractFASGAI
#' 
#' @author Xiao Nan <\url{http://www.road2stat.com}>
#' 
#' @export extractFASGAI
#' 
#' @references
#' TBD
#' 
#' @examples
#' x = readFASTA(system.file('protseq/P00750.fasta', package = 'protr'))[[1]]
#' # extractFASGAI(x)
#' 

extractFASGAI = function (x, propmat, factors, scores, lag, silent) {
  
  if (protcheck(x) == FALSE) stop('x has unrecognized amino acid type')
  
  factors = min(factors, ncol(propmat), nrow(propmat))
  
  prop.fa = factanal(propmat, factors = factors, scores = scores)
  prop.pred = predict(prop.fa)
  
  return(NULL)
  
}

