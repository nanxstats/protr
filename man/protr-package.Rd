\docType{package}
\name{protr-package}
\alias{protr}
\alias{protr-package}
\title{Protein Sequence Feature Extraction with R}
\description{
The protr package aims for protein sequence feature
extraction, which could be easily applied in
Chemoinformatics, Bioinforamtics and Chemogenomics
research. The features listed in this package include Amino
Acid Composition (Amino Acid Composition/Dipeptide
Composition/Tripeptide Composition), Autocorrelation
(Normalized Moreau-Broto Autocorrelation/Moran
Autocorrelation/Geary Autocorrelation), CTD
(Composition/Transition/Distribution), Conjoint Traid,
Quasi-sequence Order (Sequence Order Coupling
Number/Quasi-sequence Order Descriptors), and Pseudo Amino
Acid Composition (Pseudo Amino Acid Composition/Amphiphilic
Pseudo Amino Acid Composition). The package is developed by
Computational Biology and Drug Design Group, Central South
University.
}
\details{
\tabular{ll}{ Package: \tab protr\cr Type: \tab Package\cr
Version: \tab 0.2-0\cr License: \tab BSD 3-clause
License\cr }
}
\note{
The comprehensive user's guide could be opened with
\code{vignette('protr')}, which explains each descriptor
included in this package and corresponding usage.

The web interface for this package, \code{ProtrWeb} is
located at: \url{http://cbdd.csu.edu.cn:8080/protrweb/}.

Bug reports and feature requests should be sent to
\url{https://github.com/road2stat/protr/issues}.
}
\examples{
NULL
}
\author{
Nan Xiao <\email{road2stat@gmail.com}> Qingsong Xu
<\email{dasongxu@gmail.com}> Dongsheng Cao
<\email{oriental-cds@163.com}>
}
\references{
Nan Xiao, Qingsong Xu and Dongsheng Cao. protr: a
comprehensive, flexible, open-source R package for protein
sequence feature extraction. (to appear)
}
\seealso{
See \code{\link{extractAAC}}, \code{\link{extractDC}},
\code{\link{extractTC}}, \code{\link{extractMoreauBroto}},
\code{\link{extractMoran}}, \code{\link{extractGeary}},
\code{\link{extractCTDC}}, \code{\link{extractCTDT}},
\code{\link{extractCTDD}}, \code{\link{extractCTriad}},
\code{\link{extractSOCN}}, \code{\link{extractQSO}},
\code{\link{extractPAAC}}, \code{\link{extractAPAAC}} for
details.
}
\keyword{acid}
\keyword{amino}
\keyword{bioinforamtics}
\keyword{chemogenomics}
\keyword{chemoinformatics}
\keyword{descriptors}
\keyword{extraction}
\keyword{feature}
\keyword{package}
\keyword{protein}
\keyword{protr}
\keyword{sequence}

