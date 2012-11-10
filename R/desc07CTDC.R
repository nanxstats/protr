#' CTD Descriptors - Composition
#'
#' CTD Descriptors - Composition
#' 
#' This function calculates the Composition descriptor of the 
#' CTD descriptors (Dim: 21).
#' 
#' @param x A character vector, as the input protein sequence. 
#'
#' @return A length 21 named vector
#' 
#' @keywords extract CTD CTDC extractCTDC Composition
#'
#' @aliases extractCTDC
#' 
#' @author Xiao Nan <\url{http://www.road2stat.com}>
#' 
#' @seealso See \code{\link{extractCTDT}} and \code{\link{extractCTDD}} 
#'          for Transition and Distribution of the CTD descriptors.
#' 
#' @export extractCTDC
#' 
#' @references
#' Inna Dubchak, Ilya Muchink, Stephen R. Holbrook and Sung-Hou Kim. 
#' Prediction of protein folding class using global description of 
#' amino acid sequence. \emph{Proceedings of the National Academy of Sciences}. 
#' USA, 1995, 92, 8700-8704.
#' 
#' Inna Dubchak, Ilya Muchink, Christopher Mayor, Igor Dralyuk and Sung-Hou Kim.
#' Recognition of a Protein Fold in the Context of the SCOP classification. 
#' \emph{Proteins: Structure, Function and Genetics}, 1999, 35, 401-407.
#' 
#' @examples
#' x = readFASTA(system.file('AAseq/P00750.fasta', package = 'rdpi'))[[1]]
#' extractCTDC(x)
#' 

extractCTDC = function (x) {
  
  group1 = list(hydrophobicity  = c('R', 'K', 'E', 'D', 'Q', 'N'),
                normwaalsvolume = c('G', 'A', 'S', 'T', 'P', 'D', 'C'),
                polarity        = c('L', 'I', 'F', 'W', 'C', 'M', 'V', 'Y'),
                polarizability  = c('G', 'A', 'S', 'D', 'T'),
                charge          = c('K', 'R'),
                secondarystruct = c('E', 'A', 'L', 'M', 'Q', 'K', 'R', 'H'),
                solventaccess   = c('A', 'L', 'F', 'C', 'G', 'I', 'V', 'W'))
  
  group2 = list(hydrophobicity  = c('G', 'A', 'S', 'T', 'P', 'H', 'Y'),
                normwaalsvolume = c('N', 'V', 'E', 'Q', 'I', 'L'),
                polarity        = c('P', 'A', 'T', 'G', 'S'),
                polarizability  = c('C', 'P', 'N', 'V', 'E', 'Q', 'I', 'L'),
                charge          = c('A', 'N', 'C', 'Q', 'G', 'H', 'I', 'L', 
                                    'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'),
                secondarystruct = c('V', 'I', 'Y', 'C', 'W', 'F', 'T'),
                solventaccess   = c('R', 'K', 'Q', 'E', 'N', 'D'))
  
  group3 = list(hydrophobicity  = c('C', 'L', 'V', 'I', 'M', 'F', 'W'),
                normwaalsvolume = c('M', 'H', 'K', 'F', 'R', 'Y', 'W'),
                polarity        = c('H', 'Q', 'R', 'K', 'N', 'E', 'D'),
                polarizability  = c('K', 'M', 'H', 'F', 'R', 'Y', 'W'),
                charge          = c('D', 'E'),
                secondarystruct = c('G', 'N', 'P', 'S', 'D'),
                solventaccess   = c('M', 'S', 'P', 'T', 'H', 'Y'))
  
  xSplitted = strsplit(x, split = '')[[1]]
  n  = nchar(x)
  
  G = vector('list', 7)
  for (i in 1:7) G[[i]] = rep(NA, n)
  
  # Get groups for each property & each amino acid
  
  for (i in 1:7) {
    try(G[[i]][which(xSplitted %in% group1[[i]])] <- 'G1')
    try(G[[i]][which(xSplitted %in% group2[[i]])] <- 'G2')
    try(G[[i]][which(xSplitted %in% group3[[i]])] <- 'G3')
  }
  
  G = lapply(G, as.factor)
  
  CTDC = unlist(lapply(G, summary))/n
  names(CTDC) = paste('prop', rep(1:7, each = 3), 
                      '.', names(CTDC), sep = '')
  
  return(CTDC)
  
}

