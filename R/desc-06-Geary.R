#' Geary Autocorrelation Descriptor
#'
#' Geary Autocorrelation Descriptor
#'
#' This function calculates the Geary
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
#' @param customprops A \code{n x 21} named data frame contains \code{n}
#'                    customize property. Each row contains one property.
#'                    The column order for different amino acid types is
#'                    \code{'AccNo'}, \code{'A'}, \code{'R'}, \code{'N'},
#'                    \code{'D'}, \code{'C'}, \code{'E'}, \code{'Q'},
#'                    \code{'G'}, \code{'H'}, \code{'I'}, \code{'L'},
#'                    \code{'K'}, \code{'M'}, \code{'F'}, \code{'P'},
#'                    \code{'S'}, \code{'T'}, \code{'W'}, \code{'Y'},
#'                    \code{'V'}, and the columns should also be \emph{exactly}
#'                    named like this.
#'                    The \code{AccNo} column contains the properties' names.
#'                    Then users should explicitly specify these properties
#'                    with these names in the argument \code{props}.
#'                    See the examples below for a demonstration.
#'                    The default value for \code{customprops} is \code{NULL}.
#'
#' @return A length \code{nlag} named vector
#'
#' @keywords extract Geary autocorrelation
#'
#' @aliases extractGeary
#'
#' @author Nan Xiao <\url{http://nanx.me}>
#'
#' @seealso See \code{\link{extractMoreauBroto}} and \code{\link{extractMoran}}
#'          for Moreau-Broto autocorrelation descriptors and
#'          Moran autocorrelation descriptors.
#'
#' @export extractGeary
#'
#' @note For this descriptor type, users need to intelligently evaluate
#' the underlying details of the descriptors provided, instead of using
#' this function with their data blindly. It would be wise to use some
#' negative and positive control comparisons where relevant to help guide
#' interpretation of the results.
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
#' x = readFASTA(system.file('protseq/P00750.fasta', package = 'protr'))[[1]]
#' extractGeary(x)
#'
#' myprops = data.frame(AccNo = c("MyProp1", "MyProp2", "MyProp3"),
#'                      A = c(0.62,  -0.5, 15),  R = c(-2.53,   3, 101),
#'                      N = c(-0.78,  0.2, 58),  D = c(-0.9,    3, 59),
#'                      C = c(0.29,    -1, 47),  E = c(-0.74,   3, 73),
#'                      Q = c(-0.85,  0.2, 72),  G = c(0.48,    0, 1),
#'                      H = c(-0.4,  -0.5, 82),  I = c(1.38, -1.8, 57),
#'                      L = c(1.06,  -1.8, 57),  K = c(-1.5,    3, 73),
#'                      M = c(0.64,  -1.3, 75),  F = c(1.19, -2.5, 91),
#'                      P = c(0.12,     0, 42),  S = c(-0.18, 0.3, 31),
#'                      T = c(-0.05, -0.4, 45),  W = c(0.81, -3.4, 130),
#'                      Y = c(0.26,  -2.3, 107), V = c(1.08, -1.5, 43))
#'
#' # Use 4 properties in the AAindex database, and 3 cutomized properties
#' extractGeary(x, customprops = myprops,
#'              props = c('CIDH920105', 'BHAR880101',
#'                        'CHAM820101', 'CHAM820102',
#'                        'MyProp1', 'MyProp2', 'MyProp3'))

extractGeary = function (x, props = c('CIDH920105', 'BHAR880101',
                                      'CHAM820101', 'CHAM820102',
                                      'CHOC760101', 'BIGC670101',
                                      'CHAM810101', 'DAYM780201'),
                         nlag = 30L, customprops = NULL) {

  if (protcheck(x) == FALSE) stop('x has unrecognized amino acid type')

  if (nchar(x) <= nlag) {
    warning('extractGeary: length of the sequence is <= nlag; NAs will be generated')
  }

  # 1. Compute Pr values for each type of property

  AAidx = read.csv(system.file('sysdata/AAidx.csv', package = 'protr'), header = TRUE)

  if (!is.null(customprops)) AAidx = rbind(AAidx, customprops)

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

  # 3. Compute Geary Autocorrelation Descriptor

  Geary = vector('list', n)
  N  = length(xSplitted)

  Pbar = sapply(P, mean)

  for (i in 1:n) {
    for (j in 1:nlag) {
      Geary[[i]][j] = ifelse(N - j > 0,
                             ((N - 1)/(2 * (N - j))) * ((sum((P[[i]][1:(N - j)] - P[[i]][(1:(N - j)) + j])^2))/(sum((P[[i]] - Pbar[i])^2))),
                             NA)
    }
  }

  Geary = unlist(Geary)

  names(Geary) = as.vector(t(outer(props,
                                   paste('.lag', 1:nlag, sep = ''),
                                   paste, sep = '')))

  return(Geary)

}
