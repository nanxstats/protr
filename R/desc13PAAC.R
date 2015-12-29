#' Pseudo Amino Acid Composition Descriptor
#'
#' Pseudo Amino Acid Composition Descriptor
#'
#' This function calculates the Pseudo Amino Acid Composition (PAAC) descriptor
#' (Dim: \code{20 + lambda}, default is 50).
#'
#' @param x A character vector, as the input protein sequence.
#'
#' @param props A character vector, specifying the properties used.
#'              3 properties are used by default, as listed below:
#'              \describe{
#'              \item{\code{'Hydrophobicity'}}{Hydrophobicity value
#'              of the 20 amino acids}
#'              \item{\code{'Hydrophilicity'}}{Hydrophilicity value
#'              of the 20 amino acids}
#'              \item{\code{'SideChainMass'}}{Side-chain mass
#'              of the 20 amino acids}}
#'
#' @param lambda The lambda parameter for the PAAC descriptors, default is 30.
#'
#' @param w The weighting factor, default is 0.05.
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
#' @return A length \code{20 + lambda} named vector
#'
#' @note Note the default \code{20 * 3} \code{prop} values have been already
#'       independently given in the function. Users could also specify
#'       other (up to 544) properties with the Accession Number in
#'       the \code{\link{AAindex}} data, with or without the default
#'       three properties, which means users should explicitly specify
#'       the properties to use. For this descriptor type, users need to
#'       intelligently evaluate the underlying details of the descriptors
#'       provided, instead of using this function with their data blindly.
#'       It would be wise to use some negative and positive control comparisons
#'       where relevant to help guide interpretation of the results.
#'
#' @keywords extract PAAC Pseudo Composition
#'
#' @aliases extractPAAC
#'
#' @author Nan Xiao <\url{http://nanx.me}>
#'
#' @seealso See \code{\link{extractAPAAC}} for amphiphilic pseudo
#'          amino acid composition descriptor.
#'
#' @export extractPAAC
#'
#' @references
#' Kuo-Chen Chou. Prediction of Protein Cellular Attributes
#' Using Pseudo-Amino Acid Composition.
#' \emph{PROTEINS: Structure, Function, and Genetics}, 2001, 43: 246-255.
#'
#' Kuo-Chen Chou. Using Amphiphilic Pseudo Amino Acid Composition
#' to Predict Enzyme Subfamily Classes. \emph{Bioinformatics}, 2005, 21, 10-19.
#'
#' JACS, 1962, 84: 4240-4246. (C. Tanford). (The hydrophobicity data)
#'
#' PNAS, 1981, 78:3824-3828 (T.P.Hopp & K.R.Woods). (The hydrophilicity data)
#'
#' CRC Handbook of Chemistry and Physics, 66th ed.,
#' CRC Press, Boca Raton, Florida (1985). (The side-chain mass data)
#'
#' R.M.C. Dawson, D.C. Elliott, W.H. Elliott, K.M. Jones,
#' Data for Biochemical Research 3rd ed., Clarendon Press Oxford (1986).
#' (The side-chain mass data)
#'
#' @examples
#' x = readFASTA(system.file('protseq/P00750.fasta', package = 'protr'))[[1]]
#' extractPAAC(x)
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
#' # Use 3 default properties, 4 properties in the AAindex database,
#' # and 3 cutomized properties
#' extractPAAC(x, customprops = myprops,
#'             props = c('Hydrophobicity', 'Hydrophilicity', 'SideChainMass',
#'                       'CIDH920105', 'BHAR880101',
#'                       'CHAM820101', 'CHAM820102',
#'                       'MyProp1', 'MyProp2', 'MyProp3'))

extractPAAC = function (x, props = c('Hydrophobicity', 'Hydrophilicity', 'SideChainMass'),
                        lambda = 30, w = 0.05, customprops = NULL) {

  if (protcheck(x) == FALSE) stop('x has unrecognized amino acid type')

  if (nchar(x) <= lambda) stop('Length of the protein sequence must be greater than "lambda"')

  AAidx = read.csv(system.file('sysdata/AAidx.csv', package = 'protr'), header = TRUE)

  tmp = data.frame(AccNo = c("Hydrophobicity", "Hydrophilicity", "SideChainMass"),
                   A = c(0.62,  -0.5, 15),  R = c(-2.53,   3, 101),
                   N = c(-0.78,  0.2, 58),  D = c(-0.9,    3, 59),
                   C = c(0.29,    -1, 47),  E = c(-0.74,   3, 73),
                   Q = c(-0.85,  0.2, 72),  G = c(0.48,    0, 1),
                   H = c(-0.4,  -0.5, 82),  I = c(1.38, -1.8, 57),
                   L = c(1.06,  -1.8, 57),  K = c(-1.5,    3, 73),
                   M = c(0.64,  -1.3, 75),  F = c(1.19, -2.5, 91),
                   P = c(0.12,     0, 42),  S = c(-0.18, 0.3, 31),
                   T = c(-0.05, -0.4, 45),  W = c(0.81, -3.4, 130),
                   Y = c(0.26,  -2.3, 107), V = c(1.08, -1.5, 43))

  AAidx = rbind(AAidx, tmp)

  if (!is.null(customprops)) AAidx = rbind(AAidx, customprops)

  aaidx = AAidx[, -1]
  row.names(aaidx) = AAidx[, 1]

  n = length(props)

  # Standardize H0 to H

  H0 = as.matrix(aaidx[props, ])

  H  = matrix(ncol = 20, nrow = n)
  for (i in 1:n) H[i, ] = (H0[i, ] - mean(H0[i, ]))/(sqrt(sum((H0[i, ] - mean(H0[i, ]))^2)/20))
  AADict = c('A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I',
             'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V')
  dimnames(H) = list(props, AADict)

  # Compute (big) Theta

  Theta = vector('list', lambda)

  xSplitted = strsplit(x, split = '')[[1]]

  N = length(xSplitted)

  for (i in 1:lambda) {
    for (j in 1:(N-i)) {
      Theta[[i]][j] = mean((H[, xSplitted[j]] - H[, xSplitted[j + i]])^2)
    }
  }

  # Compute (small) theta

  theta = sapply(Theta, mean)

  # Compute first 20 features

  fc  = summary(factor(xSplitted, levels = AADict), maxsum = 21)
  Xc1 = fc/(1 + (w * sum(theta)))
  names(Xc1) = paste('Xc1.', names(Xc1), sep = '')

  # Compute last lambda features

  Xc2 = (w * theta)/(1 + (w * sum(theta)))
  names(Xc2) = paste('Xc2.lambda.', 1:lambda, sep = '')

  # Combine (20 + lambda) features

  Xc = c(Xc1, Xc2)

  return(Xc)

}
