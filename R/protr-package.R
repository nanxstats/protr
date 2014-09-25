#' Protein Sequence Descriptor Calculation and Similarity Computation with R
#'
#' The protr package focus on offering a unique and comprehensive
#' toolkit for generating various numerical representation schemes of protein
#' sequence. The descriptors included in the protr package are extensively
#' utilized in Bioinformatics and Chemogenomics research. The commonly used
#' descriptors listed in protr include amino acid composition,
#' autocorrelation, CTD, conjoint traid, quasi-sequence order, pseudo amino
#' acid composition, and profile-based descriptors derived by
#' Position-Specific Scoring Matrix (PSSM). The descriptors for
#' proteochemometric (PCM) modeling, includes the generalized scales-based
#' descriptors derived by principal components analysis, factor analysis,
#' multidimensional scaling, amino acid properties (AAindex), 20+ classes of
#' 2D and 3D molecular descriptors (Topological, WHIM, VHSE, etc.), and
#' generalized BLOSUM/PAM matrix-derived descriptors. The protr package also
#' integrates the function of parallelized similarity computation derived by
#' pairwise protein sequence alignment and Gene Ontology (GO) semantic
#' similarity measures. ProtrWeb, the web server built on protr, is located
#' at: http://cbdd.csu.edu.cn:8080/protrweb/
#'
#' \tabular{ll}{ Package: \tab protr\cr Type: \tab Package\cr
#' Version: \tab 0.2-0\cr License: \tab BSD 3-clause License\cr }
#'
#' @name protr-package
#' @aliases protr
#' @docType package
#' @exportPattern "^[^\\.]"
#' @author Nan Xiao <\email{road2stat@@gmail.com}>
#'         Qingsong Xu <\email{dasongxu@@gmail.com}>
#'         Dongsheng Cao <\email{oriental-cds@@163.com}>
#'
#' @note
#' The comprehensive user's guide could be opened with \code{vignette('protr')},
#' which explains each descriptor included in this package and corresponding
#' usage.
#'
#' The web server for this package, \code{ProtrWeb} is located at:
#' \url{http://cbdd.csu.edu.cn:8080/protrweb/}.
#'
#' Bug reports and feature requests should be sent to
#' \url{https://github.com/road2stat/protr/issues}.
#'
#' @references
#' (to appear)
#'
#' @keywords protr package protein sequence amino acid feature extraction
#'           descriptors chemoinformatics bioinforamtics chemogenomics
#'
#' @example inst/examples/protr-package-ex.R
NULL
