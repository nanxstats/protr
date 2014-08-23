#' Protein Sequence Descriptor Calculation and Similarity Computation with R
#'
#' The protr package focus on offering a unique and comprehensive toolkit 
#' for protein sequence descriptor calculation and similarity computation. 
#' The descriptors included in the protr package are extensively utilized 
#' in Bioinformatics and Chemogenomics research. The qualitative descriptors 
#' listed in protr include Amino Acid Composition (Amino Acid 
#' Composition/Dipeptide Composition/Tripeptide Composition) descriptor, 
#' Autocorrelation (Normalized Moreau-Broto Autocorrelation/Moran 
#' Autocorrelation/Geary Autocorrelation) descriptor, CTD 
#' (Composition/Transition/Distribution) descriptor, Conjoint Traid descriptor, 
#' Quasi-sequence Order (Sequence Order Coupling Number/Quasi-sequence Order 
#' Descriptors) descriptor and Pseudo Amino Acid Composition (Pseudo Amino 
#' Acid Composition/Amphiphilic Pseudo Amino Acid Composition) descriptor. 
#' The quantitative descriptors, for Proteochemometric (PCM) Modeling, 
#' includes the Generalized Scales-Based Descriptors derived by Principal 
#' Components Analysis, Generalized Scales-Based Descriptors derived by 
#' AA-Properties (AAindex), Generalized Scales-Based Descriptors derived 
#' by 20+ classes of 2D and 3D Molecular Descriptors (Topological, 
#' WHIM, VHSE, etc.), Generalized Scales-Based Descriptors derived by 
#' Factor Analysis, Generalized Scales-Based Descriptors derived by 
#' Multidimensional Scaling, and Generalized BLOSUM/PAM Matrix-Derived 
#' Descriptors. The protr package also integrates the functionality of 
#' parallellized similarity computation derived by protein sequence alignment 
#' and Gene Ontology (GO) semantic similarity measures between a list of 
#' protein sequences / GO terms / Entrez Gene IDs. ProtrWeb, the web service 
#' built on protr, is located at: http://cbdd.csu.edu.cn:8080/protrweb/ . 
#' The protr package is developed by Computational Biology and Drug Design 
#' (CBDD) Group, Central South University.
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
