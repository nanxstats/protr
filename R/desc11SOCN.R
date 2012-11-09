#' Sequence-Order-Coupling Numbers
#'
#' Sequence-Order-Coupling Numbers
#' 
#' This function calculates the Sequence-Order-Coupling Numbers 
#' (Dim: \code{nlag}).
#' 
#' @param x A character vector, as the input protein sequence. 
#' 
#' @param dist A character, specify the distance matrix used between 
#'             the 20 amini acids. One of 'Schneider-Wrede' or 'Grantham'
#'
#' @param nlag The maximum lag, defualt is 60.
#'
#' @return A length \code{nlag} named vector
#' 
#' @keywords extract SOCN extractSOCN Sequence Order Coupling Number
#'
#' @aliases extractSOCN
#' 
#' @author Xiao Nan <\url{http://www.road2stat.com}>
#' 
#' @seealso See \code{\link{extractQSO}} for quasi-sequence-order descriptors.
#' 
#' @export extractSOCN
#' 
#' @references
#' Kuo-Chen Chou. Prediction of Protein Subcellar Locations by 
#' Incorporating Quasi-Sequence-Order Effect. 
#' \emph{Biochemical and Biophysical Research Communications}, 
#' 2000, 278, 477-483.
#' 
#' Gisbert Schneider and Paul Wrede. The Rational Design of 
#' Amino Acid Sequences by Artificial Neural Networks and 
#' Simulated Molecular Evolution: Do Novo Design of an Idealized 
#' Leader Cleavage Site. 
#' \emph{Biophys Journal}, 1994, 66, 335-344.
#' 
#' Grantham, R. Amino acid difference formula to help explain protein evolution.
#' \emph{Science}, 1974, 185, 862-864.
#' 
#' @examples
#' A06852 = readFASTA(system.file('AAseq/A06852.fasta', package = 'rdpi'))
#' x = 'MPRLFSYLLGVWLLLSQLPREIPGQSTNDFIKACGRELVRLWVEICGSVSWGRTALSLEEPQLETGPPAETMPSSITKDAEILKMMLEFVPNLPQELKATLSERQPSLRELQQSASKDSNLNFEEFKKIILNRQNEAEDKSLLELKNLGLDKHSRKKRLFRMTLSEKCCQVGCIRKDIARLC'
#' extractSOCN(x, dist = 'Schneider-Wrede')
#' 

extractSOCN = function (x, dist = c('Schneider-Wrede', 'Grantham'), nlag = 60) {
  
  dist = match.arg(dist)
  
  if (dist == 'Schneider-Wrede') {
    DistMat = read.csv(system.file('Schneider-Wrede.csv', package = 'rdpi'), 
                       header = TRUE)
  } else if (dist == 'Grantham') {
    DistMat = read.csv(system.file('Grantham.csv', package = 'rdpi'), 
                       header = TRUE)
  } else {
    stop("The 'dist' argument must be one of 'Schneider-Wrede' or 'Grantham'")
  }
  
  row.names(DistMat) = as.character(DistMat[, 1])
  DistMat = DistMat[, -1]
  
  xSplitted = strsplit(x, split = '')[[1]]
  N = nchar(x)
  
  # Compute tau_d
  
  tau = vector('list', nlag)
  
  for (d in 1:nlag) {
    for (i in 1:(N - d)) {
      tau[[d]][i] = (DistMat[xSplitted[i], xSplitted[i + d]])^2
    }
  }
  
  tau = sapply(tau, sum)
  
  names(tau) = paste('tau', 1:nlag, sep = '')
  
}

