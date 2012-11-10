#' Moran Autocorrelation Descriptor
#'
#' Moran Autocorrelation Descriptor
#' 
#' This function calculates the Moran
#' autocorrelation descriptor (Dim: \code{length(props) * nlag}).
#' 
#' @param x A character vector, as the input protein sequence. 
#' 
#' @param props A character vector, specifying the 
#'              Accession Number of the target properties. 
#'              8 properties are used by default, as listed below:
#'              \describe{
#'              \item{AccNo. CIDH920105}{Normalized average hydrophobicity scales (Cid et al., 1992)}
#'              \item{AccNo. BHAR880101}{Average flexibility indices (Bhaskaran-Ponnuswamy, 1988)}
#'              \item{AccNo. CHAM820101}{Polarizability parameter (Charton-Charton, 1982)}
#'              \item{AccNo. CHAM820102}{Free energy of solution in water, kcal/mole (Charton-Charton, 1982)}
#'              \item{AccNo. CHOC760101}{Residue accessible surface area in tripeptide (Chothia, 1976)}
#'              \item{AccNo. BIGC670101}{Residue volume (Bigelow, 1967)}
#'              \item{AccNo. CHAM810101}{Steric parameter (Charton, 1981)}
#'              \item{AccNo. DAYM780201}{Relative mutability (Dayhoff et al., 1978b)}}
#'              
#' @param nlag Maximum value of the lag parameter. Default is \code{30}.
#'
#' @return A length \code{nlag} named vector
#' 
#' @keywords extract extractMoran Moran autocorrelation
#'
#' @aliases extractMoran
#' 
#' @author Xiao Nan <\url{http://www.road2stat.com}>
#' 
#' @seealso See \code{\link{extractMoreauBroto}} and \code{\link{extractGeary}} 
#'          for Moreau-Broto autocorrelation descriptors and 
#'          Geary autocorrelation descriptors.
#' 
#' @export extractMoran
#' 
#' @references
#' AAindex: Amino acid index database.
#' \url{http://www.genome.ad.jp/dbget/aaindex.html}
#' 
#' Feng, Z.P. and Zhang, C.T. (2000)
#' Prediction of membrane protein types based on the hydrophobic
#' index of amino acids.
#' \emph{Journal of Protein Chemistry}, 19, 269-275.
#' 
#' Horne, D.S. (1988)
#' Prediction of protein helix content from
#' an autocorrelation analysis of sequence hydrophobicities.
#' \emph{Biopolymers}, 27, 451-477.
#' 
#' Sokal, R.R. and Thomson, B.A. (2006)
#' Population structure inferred by local spatial autocorrelation: 
#' an Usage from an Amerindian tribal population. 
#' \emph{American Journal of Physical Anthropology}, 129, 121-131.
#' 
#' @examples
#' x = readFASTA(system.file('AAseq/P00750.fasta', package = 'rdpi'))[[1]]
#' extractMoran(x)
#' 

extractMoran = function (x, props = c('CIDH920105', 'BHAR880101',
                                      'CHAM820101', 'CHAM820102',
                                      'CHOC760101', 'BIGC670101',
                                      'CHAM810101', 'DAYM780201'), 
                         nlag = 30L) {
  
  # 1. Compute Pr values for each type of property
  
  AAidx = read.csv(system.file('AAidx.csv', package = 'rdpi'), header = TRUE)
  aaidx = AAidx[, -1]
  row.names(aaidx) = AAidx[, 1]
  n = length(props)
  pmean = rowMeans(aaidx[props, ])
  psd   = apply(aaidx[props, ], 1, sd) * sqrt((20 - 1)/20) # sd() uses (n-1)
  
  Pr = data.frame(matrix(ncol = 20, nrow = n))
  for (i in 1:n) Pr[i, ] = (aaidx[props[i], ] - pmean[i])/psd[i]
  
  # 2. Replace character with numbers, also applies to less than 20 AA occured
  
  xSplitted = strsplit(x, split = '')[[1]]
  AADict = c('A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 
             'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V')
  names(Pr) = AADict
  
  P = vector('list', n)
  for (i in 1:n) P[[i]] = xSplitted
  
  for (i in 1:n) {
    for (j in AADict) {
      try(P[[i]][which(P[[i]] == j)] <- Pr[i, j], silent = TRUE)
    }
  }
  
  P = lapply(P, as.numeric)
  
  # 3. Compute Moran Autocorrelation Descriptor
  
  Moran = vector('list', n)
  N  = length(xSplitted)
  
  Pbar = sapply(P, mean)
  
  for (i in 1:n) {
    for (j in 1:nlag) {
      Moran[[i]][j] = (N/(N - j)) * ((sum((P[[i]][1:(N - j)] - Pbar[i]) * (P[[i]][(1:(N - j)) + j] - Pbar[i])))/(sum((P[[i]] - Pbar[i])^2)))
    }
  }
  
  Moran = unlist(Moran)
  
  names(Moran) = as.vector(t(outer(props, 
                                   paste('.lag', 1:nlag, sep = ''), 
                                   paste, sep = '')))
  
  return(Moran)
  
}

