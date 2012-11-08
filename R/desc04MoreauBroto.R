#' Normalized Moreau-Broto Autocorrelation Descriptor
#'
#' Normalized Moreau-Broto Autocorrelation Descriptor
#' 
#' This function calculates the normalized Moreau-Broto 
#' autocorrelation descriptor (total \code{length(props) * nlag}).
#' 
#' @param x A character vector, as the input protein sequence. 
#' 
#' @param props A character vector, specifying the 
#'              Accession Number of the target properties. 
#'              8 properties are used by default.
#'              
#' @param nlag Maximum value of the lag parameter. Default is \code{30}.
#'
#' @return A length \code{nlag} named vector
#' 
#' @keywords extract extractMoreauBroto normalized autocorrelation
#'           MoreauBroto Moreau-Broto Moreau Broto
#'
#' @aliases extractMoreauBroto
#' 
#' @author Xiao Nan <\url{http://www.road2stat.com}>
#' 
#' @seealso See \code{\link{extractMoran}} and \code{\link{extractGeary}} 
#'          for Moran autocorrelation descriptors and 
#'          Geary autocorrelation descriptors.
#' 
#' @export extractMoreauBroto
#' 
#' @references
#' References needed here!
#' 
#' @examples
#' A06852 = readFASTA(system.file('AAseq/A06852.fasta', package = 'rdpi'))
#' x = 'MPRLFSYLLGVWLLLSQLPREIPGQSTNDFIKACGRELVRLWVEICGSVSWGRTALSLEEPQLETGPPAETMPSSITKDAEILKMMLEFVPNLPQELKATLSERQPSLRELQQSASKDSNLNFEEFKKIILNRQNEAEDKSLLELKNLGLDKHSRKKRLFRMTLSEKCCQVGCIRKDIARLC'
#' extractMoreauBroto(x)
#' 

extractMoreauBroto = function (x, props = c('ANDN920101', 'ARGP820101', 
                                            'ARGP820102', 'ARGP820103', 
                                            'BEGF750101', 'BEGF750102', 
                                            'BEGF750103', 'BHAR880101'), 
                               nlag = 30L) {
  
  # 1. Compute Pr values for each type of property
  
  AAidx = read.csv(system.file('AAidx.csv', package = 'rdpi'), header = TRUE)
  aaidx = AAidx[, -1]
  row.names(aaidx) = AAidx[, 1]
  n = length(props)
  # is this *n* for sigma right?
  pmean = rowMeans(aaidx[props, ])
  psd   = apply(aaidx[props, ], 1, sd) * sqrt((20 - 1)/20)
  
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
  
  # 3. Compute Moreau-Broto Descriptor
  
  MB = vector('list', n)
  N  = length(xSplitted)
  
  for (i in 1:n) {
    for (j in 1:nlag) {
      MB[[i]][j] = sum(P[[i]][1:(N - j)] * P[[i]][(1:(N - j)) + j])/(N - j)
    }
  }
  
  MB = unlist(MB)
  
  names(MB) = as.vector(t(outer(props, 
                                paste('Lag', 1:nlag, sep = ''), 
                                paste, sep = '')))
  
  return(MB)
  
}

