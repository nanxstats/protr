#' Generating Various Numerical Representation Schemes of Protein Sequence
#'
#' The protr package is a comprehensive toolkit for
#' generating various numerical representation schemes of protein sequence.
#' The descriptors are extensively utilized in bioinformatics and
#' chemogenomics research. The commonly used descriptors include amino acid
#' composition, autocorrelation, CTD, conjoint traid, quasi-sequence order,
#' pseudo amino acid composition, and profile-based descriptors derived by
#' Position-Specific Scoring Matrix (PSSM). The descriptors for
#' proteochemometric (PCM) modeling include the scales-based
#' descriptors derived by principal components analysis, factor analysis,
#' multidimensional scaling, amino acid properties (AAindex), 20+ classes of
#' 2D and 3D molecular descriptors (Topological, WHIM, VHSE, etc.), and
#' BLOSUM/PAM matrix-derived descriptors. The protr package also
#' integrates the function of parallelized similarity computation derived by
#' pairwise protein sequence alignment and Gene Ontology (GO) semantic
#' similarity measures.
#'
#' \tabular{ll}{ Package: \tab protr\cr Type: \tab Package\cr
#' Version: \tab 1.1-0\cr License: \tab BSD 3-clause License\cr }
#'
#' @name protr-package
#' @aliases protr
#' @docType package
#' @author Nan Xiao <\email{road2stat@@gmail.com}>
#'         Qing-Song Xu <\email{dasongxu@@gmail.com}>
#'         Dong-Sheng Cao <\email{oriental-cds@@163.com}>
#'
#' @note
#' The user guide can be opened with \code{vignette('protr')},
#' which explains every descriptor included in this package and
#' their corresponding usage.
#'
#' The web server for this package, \code{ProtrWeb} is located at:
#' \url{http://protr.org}.
#'
#' Bug reports and feature requests should be sent to
#' \url{https://github.com/road2stat/protr/issues}.
#'
#' @references
#' Xiao, N., Cao, D.-S., Zhu, M.-F., and Xu, Q.-S. (2015).
#' protr/ProtrWeb: R package and web server for generating various
#' numerical representation schemes of protein sequences.
#' \emph{Bioinformatics} 31, 1857--1859.
#'
#' @importFrom stats cmdscale dist factanal na.omit prcomp predict sd
#' @importFrom utils combn read.csv
#'
#' @example inst/examples/protr-package-ex.R
NULL
