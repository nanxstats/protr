#' AAindex Data of 544 Physicochemical and Biological Properties
#' for 20 Amino Acids
#'
#' The data was extracted from the AAindex1 database ver 9.1
#' (\url{ftp://ftp.genome.jp/pub/db/community/aaindex/aaindex1}) as of
#' November, 2012 (Data Last Modified 2006-08-14).
#'
#' With this dataset, users can investigate each property's accession number
#' and other details. Visit \url{http://www.genome.jp/dbget/aaindex.html}
#' for more information.
#'
#' @docType data
#' @name AAindex
#' @aliases AAindex
#' @keywords AAindex datasets
#' @examples
#' data(AAindex)
#'
NULL

#' OptAA3d.sdf - 20 Amino Acids Optimized with MOE 2011.10 (Semiempirical AM1)
#'
#' OptAA3d.sdf - 20 Amino Acids Optimized with MOE 2011.10 (Semiempirical AM1)
#'
#' @docType data
#' @name OptAA3d
#' @aliases OptAA3d
#' @keywords datasets
#' @examples
#' # This operation requires the rcdk package
#' # require(rcdk)
#' # optaa3d = load.molecules(system.file("sysdata/OptAA3d.sdf", package = "protr"))
#' # view.molecule.2d(optaa3d[[1]])  # view the first AA
NULL

#' Meta Information for the 20 Amino Acids
#'
#' This dataset includes the meta information of
#' the 20 amino acids used for the 2D and 3D descriptor
#' calculation in this package. Each column represents:
#' \itemize{
#' \item{\code{AAName}} {Amino Acid Name}
#' \item{\code{Short}} {One-Letter Representation}
#' \item{\code{Abbreviation}} {Three-Letter Representation}
#' \item{\code{mol}} {SMILE Representation}
#' \item{\code{PUBCHEM_COMPOUND_CID}} {PubChem CID for the Amino Acid}
#' \item{\code{PUBCHEM_LINK}} {PubChem Link for the Amino Acid}
#' }
#'
#' @docType data
#' @name AAMetaInfo
#' @aliases AAMetaInfo
#' @keywords datasets
#' @examples
#' data(AAMetaInfo)
#'
NULL

#' 2D Descriptors for 20 Amino Acids calculated by MOE 2011.10
#'
#' This dataset includes the 2D descriptors of
#' the 20 amino acids calculated by MOE 2011.10
#' used for scales extraction in this package.
#'
#' @docType data
#' @name AAMOE2D
#' @aliases AAMOE2D
#' @keywords datasets
#' @examples
#' data(AAMOE2D)
#'
NULL

#' 3D Descriptors for 20 Amino Acids calculated by MOE 2011.10
#'
#' This dataset includes the 3D descriptors of the 20 amino acids
#' calculated by MOE 2011.10 used for scales extraction in this package.
#' All amino acid molecules had also been optimized with MOE (semiempirical AM1)
#' before calculating these 3D descriptors.
#' The SDF file containing the information of the optimized amino acid molecules
#' is included in this package. See \code{\link{OptAA3d}} for more information.
#'
#' @docType data
#' @name AAMOE3D
#' @aliases AAMOE3D
#' @keywords datasets
#' @examples
#' data(AAMOE3D)
#'
NULL

#' CPSA Descriptors for 20 Amino Acids calculated by Discovery Studio
#'
#' This dataset includes the CPSA descriptors of the 20 amino acids
#' calculated by Discovery Studio (version 2.5) used for scales extraction
#' in this package.
#'
#' All amino acid molecules had also been optimized with MOE 2011.10
#' (semiempirical AM1) before calculating these CPSA descriptors.
#' The SDF file containing the information of the optimized amino acid
#' molecules is included in this package. See \code{\link{OptAA3d}}
#' for more information.
#'
#' @docType data
#' @name AACPSA
#' @aliases AACPSA
#' @keywords datasets
#' @examples
#' data(AACPSA)
#'
NULL

#' All 2D Descriptors for 20 Amino Acids calculated by Dragon
#'
#' This dataset includes all the 2D descriptors of
#' the 20 amino acids calculated by Dragon (version 5.4)
#' used for scales extraction in this package.
#'
#' @docType data
#' @name AADescAll
#' @aliases AADescAll
#' @keywords datasets
#' @examples
#' data(AADescAll)
#'
NULL



#' 2D Autocorrelations Descriptors for 20 Amino Acids calculated by Dragon
#'
#' This dataset includes the 2D autocorrelations descriptors of
#' the 20 amino acids calculated by Dragon (version 5.4)
#' used for scales extraction in this package.
#'
#' @docType data
#' @name AA2DACOR
#' @aliases AA2DACOR
#' @keywords datasets
#' @examples
#' data(AA2DACOR)
#'
NULL

#' 3D-MoRSE Descriptors for 20 Amino Acids calculated by Dragon
#'
#' This dataset includes the 3D-MoRSE descriptors of
#' the 20 amino acids calculated by Dragon (version 5.4)
#' used for scales extraction in this package.
#'
#' @docType data
#' @name AA3DMoRSE
#' @aliases AA3DMoRSE
#' @keywords datasets
#' @examples
#' data(AA3DMoRSE)
#'
NULL

#' Atom-Centred Fragments Descriptors for 20 Amino Acids calculated by Dragon
#'
#' This dataset includes the atom-centred fragments descriptors of
#' the 20 amino acids calculated by Dragon (version 5.4)
#' used for scales extraction in this package.
#'
#' @docType data
#' @name AAACF
#' @aliases AAACF
#' @keywords datasets
#' @examples
#' data(AAACF)
#'
NULL

#' Burden Eigenvalues Descriptors for 20 Amino Acids calculated by Dragon
#'
#' This dataset includes the Burden eigenvalues descriptors of
#' the 20 amino acids calculated by Dragon (version 5.4)
#' used for scales extraction in this package.
#'
#' @docType data
#' @name AABurden
#' @aliases AABurden
#' @keywords datasets
#' @examples
#' data(AABurden)
#'
NULL

#' Connectivity Indices Descriptors for 20 Amino Acids calculated by Dragon
#'
#' This dataset includes the connectivity indices descriptors of
#' the 20 amino acids calculated by Dragon (version 5.4)
#' used for scales extraction in this package.
#'
#' @docType data
#' @name AAConn
#' @aliases AAConn
#' @keywords datasets
#' @examples
#' data(AAConn)
#'
NULL

#' Constitutional Descriptors for 20 Amino Acids calculated by Dragon
#'
#' This dataset includes the constitutional descriptors of
#' the 20 amino acids calculated by Dragon (version 5.4)
#' used for scales extraction in this package.
#'
#' @docType data
#' @name AAConst
#' @aliases AAConst
#' @keywords datasets
#' @examples
#' data(AAConst)
#'
NULL

#' Edge Adjacency Indices Descriptors for 20 Amino Acids calculated by Dragon
#'
#' This dataset includes the edge adjacency indices descriptors of
#' the 20 amino acids calculated by Dragon (version 5.4)
#' used for scales extraction in this package.
#'
#' @docType data
#' @name AAEdgeAdj
#' @aliases AAEdgeAdj
#' @keywords datasets
#' @examples
#' data(AAEdgeAdj)
#'
NULL

#' Eigenvalue-Based Indices Descriptors for 20 Amino Acids calculated by Dragon
#'
#' This dataset includes the eigenvalue-based indices descriptors of
#' the 20 amino acids calculated by Dragon (version 5.4)
#' used for scales extraction in this package.
#'
#' @docType data
#' @name AAEigIdx
#' @aliases AAEigIdx
#' @keywords datasets
#' @examples
#' data(AAEigIdx)
#'
NULL

#' Functional Group Counts Descriptors for 20 Amino Acids calculated by Dragon
#'
#' This dataset includes the functional group counts descriptors of
#' the 20 amino acids calculated by Dragon (version 5.4)
#' used for scales extraction in this package.
#'
#' @docType data
#' @name AAFGC
#' @aliases AAFGC
#' @keywords datasets
#' @examples
#' data(AAFGC)
#'
NULL

#' Geometrical Descriptors for 20 Amino Acids calculated by Dragon
#'
#' This dataset includes the geometrical descriptors of
#' the 20 amino acids calculated by Dragon (version 5.4)
#' used for scales extraction in this package.
#'
#' @docType data
#' @name AAGeom
#' @aliases AAGeom
#' @keywords datasets
#' @examples
#' data(AAGeom)
#'
NULL

#' GETAWAY Descriptors for 20 Amino Acids calculated by Dragon
#'
#' This dataset includes the GETAWAY descriptors of
#' the 20 amino acids calculated by Dragon (version 5.4)
#' used for scales extraction in this package.
#'
#' @docType data
#' @name AAGETAWAY
#' @aliases AAGETAWAY
#' @keywords datasets
#' @examples
#' data(AAGETAWAY)
#'
NULL

#' Information Indices Descriptors for 20 Amino Acids calculated by Dragon
#'
#' This dataset includes the information indices descriptors of
#' the 20 amino acids calculated by Dragon (version 5.4)
#' used for scales extraction in this package.
#'
#' @docType data
#' @name AAInfo
#' @aliases AAInfo
#' @keywords datasets
#' @examples
#' data(AAInfo)
#'
NULL

#' Molecular Properties Descriptors for 20 Amino Acids calculated by Dragon
#'
#' This dataset includes the molecular properties descriptors of
#' the 20 amino acids calculated by Dragon (version 5.4)
#' used for scales extraction in this package.
#'
#' @docType data
#' @name AAMolProp
#' @aliases AAMolProp
#' @keywords datasets
#' @examples
#' data(AAMolProp)
#'
NULL

#' Randic Molecular Profiles Descriptors for 20 Amino Acids calculated by Dragon
#'
#' This dataset includes the Randic molecular profiles descriptors of
#' the 20 amino acids calculated by Dragon (version 5.4)
#' used for scales extraction in this package.
#'
#' @docType data
#' @name AARandic
#' @aliases AARandic
#' @keywords datasets
#' @examples
#' data(AARandic)
#'
NULL

#' RDF Descriptors for 20 Amino Acids calculated by Dragon
#'
#' This dataset includes the RDF descriptors of
#' the 20 amino acids calculated by Dragon (version 5.4)
#' used for scales extraction in this package.
#'
#' @docType data
#' @name AARDF
#' @aliases AARDF
#' @keywords datasets
#' @examples
#' data(AARDF)
#'
NULL

#' Topological Descriptors for 20 Amino Acids calculated by Dragon
#'
#' This dataset includes the topological descriptors of
#' the 20 amino acids calculated by Dragon (version 5.4)
#' used for scales extraction in this package.
#'
#' @docType data
#' @name AATopo
#' @aliases AATopo
#' @keywords datasets
#' @examples
#' data(AATopo)
#'
NULL

#' Topological Charge Indices Descriptors for 20 Amino Acids calculated by Dragon
#'
#' This dataset includes the topological charge indices descriptors of
#' the 20 amino acids calculated by Dragon (version 5.4)
#' used for scales extraction in this package.
#'
#' @docType data
#' @name AATopoChg
#' @aliases AATopoChg
#' @keywords datasets
#' @examples
#' data(AATopoChg)
#'
NULL

#' Walk and Path Counts Descriptors for 20 Amino Acids calculated by Dragon
#'
#' This dataset includes the walk and path counts descriptors of
#' the 20 amino acids calculated by Dragon (version 5.4)
#' used for scales extraction in this package.
#'
#' @docType data
#' @name AAWalk
#' @aliases AAWalk
#' @keywords datasets
#' @examples
#' data(AAWalk)
#'
NULL

#' WHIM Descriptors for 20 Amino Acids calculated by Dragon
#'
#' This dataset includes the WHIM descriptors of
#' the 20 amino acids calculated by Dragon (version 5.4)
#' used for scales extraction in this package.
#'
#' @docType data
#' @name AAWHIM
#' @aliases AAWHIM
#' @keywords datasets
#' @examples
#' data(AAWHIM)
#'
NULL



#' BLOSUM45 Matrix for 20 Amino Acids
#'
#' BLOSUM45 Matrix for the 20 amino acids. The matrix was extracted from the
#' \code{Biostrings} package of Bioconductor.
#'
#' @docType data
#' @name AABLOSUM45
#' @aliases AABLOSUM45
#' @keywords datasets
#' @examples
#' data(AABLOSUM45)
#'
NULL

#' BLOSUM50 Matrix for 20 Amino Acids
#'
#' BLOSUM50 Matrix for the 20 amino acids. The matrix was extracted from the
#' \code{Biostrings} package of Bioconductor.
#'
#' @docType data
#' @name AABLOSUM50
#' @aliases AABLOSUM50
#' @keywords datasets
#' @examples
#' data(AABLOSUM50)
#'
NULL

#' BLOSUM62 Matrix for 20 Amino Acids
#'
#' BLOSUM62 Matrix for the 20 amino acids. The matrix was extracted from the
#' \code{Biostrings} package of Bioconductor.
#'
#' @docType data
#' @name AABLOSUM62
#' @aliases AABLOSUM62
#' @keywords datasets
#' @examples
#' data(AABLOSUM62)
#'
NULL

#' BLOSUM80 Matrix for 20 Amino Acids
#'
#' BLOSUM80 Matrix for the 20 amino acids. The matrix was extracted from the
#' \code{Biostrings} package of Bioconductor.
#'
#' @docType data
#' @name AABLOSUM80
#' @aliases AABLOSUM80
#' @keywords datasets
#' @examples
#' data(AABLOSUM80)
#'
NULL

#' BLOSUM100 Matrix for 20 Amino Acids
#'
#' BLOSUM100 Matrix for the 20 amino acids. The matrix was extracted from the
#' \code{Biostrings} package of Bioconductor.
#'
#' @docType data
#' @name AABLOSUM100
#' @aliases AABLOSUM100
#' @keywords datasets
#' @examples
#' data(AABLOSUM100)
#'
NULL

#' PAM30 Matrix for 20 Amino Acids
#'
#' PAM30 Matrix for the 20 amino acids. The matrix was extracted from the
#' \code{Biostrings} package of Bioconductor.
#'
#' @docType data
#' @name AAPAM30
#' @aliases AAPAM30
#' @keywords datasets
#' @examples
#' data(AAPAM30)
#'
NULL

#' PAM40 Matrix for 20 Amino Acids
#'
#' PAM40 Matrix for the 20 amino acids. The matrix was extracted from the
#' \code{Biostrings} package of Bioconductor.
#'
#' @docType data
#' @name AAPAM40
#' @aliases AAPAM40
#' @keywords datasets
#' @examples
#' data(AAPAM40)
#'
NULL

#' PAM70 Matrix for 20 Amino Acids
#'
#' PAM70 Matrix for the 20 amino acids. The matrix was extracted from the
#' \code{Biostrings} package of Bioconductor.
#'
#' @docType data
#' @name AAPAM70
#' @aliases AAPAM70
#' @keywords datasets
#' @examples
#' data(AAPAM70)
#'
NULL

#' PAM120 Matrix for 20 Amino Acids
#'
#' PAM120 Matrix for the 20 amino acids. The matrix was extracted from the
#' \code{Biostrings} package of Bioconductor.
#'
#' @docType data
#' @name AAPAM120
#' @aliases AAPAM120
#' @keywords datasets
#' @examples
#' data(AAPAM120)
#'
NULL

#' PAM250 Matrix for 20 Amino Acids
#'
#' PAM250 Matrix for the 20 amino acids. The matrix was extracted from the
#' \code{Biostrings} package of Bioconductor.
#'
#' @docType data
#' @name AAPAM250
#' @aliases AAPAM250
#' @keywords datasets
#' @examples
#' data(AAPAM250)
#'
NULL
