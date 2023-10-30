#' AAindex Data of 544 Physicochemical and Biological Properties
#' for 20 Amino Acids
#'
#' The data was extracted from the AAindex1 database ver 9.1
#' (\url{https://www.genome.jp/ftp/db/community/aaindex/old/ver9.1/aaindex1})
#' as of November, 2012 (data last modified on 2006-08-14).
#'
#' With this dataset, users can investigate each property's accession number
#' and other details. See \url{https://www.genome.jp/aaindex/}
#' for more information.
#'
#' @docType data
#' @name AAindex
#' @examples
#' data(AAindex)
NULL

#' OptAA3d.sdf - 20 Amino Acids Optimized with MOE 2011.10 (Semiempirical AM1)
#'
#' OptAA3d.sdf - 20 Amino Acids Optimized with MOE 2011.10 (Semiempirical AM1)
#'
#' @docType data
#' @name OptAA3d
#' @examples
#' # this operation requires the rcdk package
#' # require(rcdk)
#' # optaa3d = load.molecules(system.file("sysdata/OptAA3d.sdf", package = "protr"))
#' # view.molecule.2d(optaa3d[[1]])  # view the first AA
NULL

#' Meta Information for the 20 Amino Acids
#'
#' This dataset includes the meta information of
#' the 20 amino acids used for the 2D and 3D descriptor
#' calculation in this package. Each column represents:
#' \describe{
#'   \item{\code{AAName}}{Amino acid name}
#'   \item{\code{Short}}{One-letter representation}
#'   \item{\code{Abbreviation}}{Three-letter representation}
#'   \item{\code{mol}}{SMILES representation}
#'   \item{\code{PUBCHEM_COMPOUND_CID}}{PubChem CID for the amino acid}
#'   \item{\code{PUBCHEM_LINK}}{PubChem link for the amino acid}
#' }
#'
#' @docType data
#' @name AAMetaInfo
#' @examples
#' data(AAMetaInfo)
NULL

#' 2D Descriptors for 20 Amino Acids calculated by MOE 2011.10
#'
#' This dataset includes the 2D descriptors of
#' the 20 amino acids calculated by MOE 2011.10
#' used for scales extraction in this package.
#'
#' @docType data
#' @name AAMOE2D
#' @examples
#' data(AAMOE2D)
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
#' @examples
#' data(AAMOE3D)
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
#' @examples
#' data(AACPSA)
NULL

#' All 2D Descriptors for 20 Amino Acids calculated by Dragon
#'
#' This dataset includes all the 2D descriptors of
#' the 20 amino acids calculated by Dragon (version 5.4)
#' used for scales extraction in this package.
#'
#' @docType data
#' @name AADescAll
#' @examples
#' data(AADescAll)
NULL



#' 2D Autocorrelations Descriptors for 20 Amino Acids calculated by Dragon
#'
#' This dataset includes the 2D autocorrelations descriptors of
#' the 20 amino acids calculated by Dragon (version 5.4)
#' used for scales extraction in this package.
#'
#' @docType data
#' @name AA2DACOR
#' @examples
#' data(AA2DACOR)
NULL

#' 3D-MoRSE Descriptors for 20 Amino Acids calculated by Dragon
#'
#' This dataset includes the 3D-MoRSE descriptors of
#' the 20 amino acids calculated by Dragon (version 5.4)
#' used for scales extraction in this package.
#'
#' @docType data
#' @name AA3DMoRSE
#' @examples
#' data(AA3DMoRSE)
NULL

#' Atom-Centred Fragments Descriptors for 20 Amino Acids calculated by Dragon
#'
#' This dataset includes the atom-centred fragments descriptors of
#' the 20 amino acids calculated by Dragon (version 5.4)
#' used for scales extraction in this package.
#'
#' @docType data
#' @name AAACF
#' @examples
#' data(AAACF)
NULL

#' Burden Eigenvalues Descriptors for 20 Amino Acids calculated by Dragon
#'
#' This dataset includes the Burden eigenvalues descriptors of
#' the 20 amino acids calculated by Dragon (version 5.4)
#' used for scales extraction in this package.
#'
#' @docType data
#' @name AABurden
#' @examples
#' data(AABurden)
NULL

#' Connectivity Indices Descriptors for 20 Amino Acids calculated by Dragon
#'
#' This dataset includes the connectivity indices descriptors of
#' the 20 amino acids calculated by Dragon (version 5.4)
#' used for scales extraction in this package.
#'
#' @docType data
#' @name AAConn
#' @examples
#' data(AAConn)
NULL

#' Constitutional Descriptors for 20 Amino Acids calculated by Dragon
#'
#' This dataset includes the constitutional descriptors of
#' the 20 amino acids calculated by Dragon (version 5.4)
#' used for scales extraction in this package.
#'
#' @docType data
#' @name AAConst
#' @examples
#' data(AAConst)
NULL

#' Edge Adjacency Indices Descriptors for 20 Amino Acids calculated by Dragon
#'
#' This dataset includes the edge adjacency indices descriptors of
#' the 20 amino acids calculated by Dragon (version 5.4)
#' used for scales extraction in this package.
#'
#' @docType data
#' @name AAEdgeAdj
#' @examples
#' data(AAEdgeAdj)
NULL

#' Eigenvalue-Based Indices Descriptors for 20 Amino Acids calculated by Dragon
#'
#' This dataset includes the eigenvalue-based indices descriptors of
#' the 20 amino acids calculated by Dragon (version 5.4)
#' used for scales extraction in this package.
#'
#' @docType data
#' @name AAEigIdx
#' @examples
#' data(AAEigIdx)
NULL

#' Functional Group Counts Descriptors for 20 Amino Acids calculated by Dragon
#'
#' This dataset includes the functional group counts descriptors of
#' the 20 amino acids calculated by Dragon (version 5.4)
#' used for scales extraction in this package.
#'
#' @docType data
#' @name AAFGC
#' @examples
#' data(AAFGC)
NULL

#' Geometrical Descriptors for 20 Amino Acids calculated by Dragon
#'
#' This dataset includes the geometrical descriptors of
#' the 20 amino acids calculated by Dragon (version 5.4)
#' used for scales extraction in this package.
#'
#' @docType data
#' @name AAGeom
#' @examples
#' data(AAGeom)
NULL

#' GETAWAY Descriptors for 20 Amino Acids calculated by Dragon
#'
#' This dataset includes the GETAWAY descriptors of
#' the 20 amino acids calculated by Dragon (version 5.4)
#' used for scales extraction in this package.
#'
#' @docType data
#' @name AAGETAWAY
#' @examples
#' data(AAGETAWAY)
NULL

#' Information Indices Descriptors for 20 Amino Acids calculated by Dragon
#'
#' This dataset includes the information indices descriptors of
#' the 20 amino acids calculated by Dragon (version 5.4)
#' used for scales extraction in this package.
#'
#' @docType data
#' @name AAInfo
#' @examples
#' data(AAInfo)
NULL

#' Molecular Properties Descriptors for 20 Amino Acids calculated by Dragon
#'
#' This dataset includes the molecular properties descriptors of
#' the 20 amino acids calculated by Dragon (version 5.4)
#' used for scales extraction in this package.
#'
#' @docType data
#' @name AAMolProp
#' @examples
#' data(AAMolProp)
NULL

#' Randic Molecular Profiles Descriptors for 20 Amino Acids calculated by Dragon
#'
#' This dataset includes the Randic molecular profiles descriptors of
#' the 20 amino acids calculated by Dragon (version 5.4)
#' used for scales extraction in this package.
#'
#' @docType data
#' @name AARandic
#' @examples
#' data(AARandic)
NULL

#' RDF Descriptors for 20 Amino Acids calculated by Dragon
#'
#' This dataset includes the RDF descriptors of
#' the 20 amino acids calculated by Dragon (version 5.4)
#' used for scales extraction in this package.
#'
#' @docType data
#' @name AARDF
#' @examples
#' data(AARDF)
NULL

#' Topological Descriptors for 20 Amino Acids calculated by Dragon
#'
#' This dataset includes the topological descriptors of
#' the 20 amino acids calculated by Dragon (version 5.4)
#' used for scales extraction in this package.
#'
#' @docType data
#' @name AATopo
#' @examples
#' data(AATopo)
NULL

#' Topological Charge Indices Descriptors for 20 Amino Acids calculated by Dragon
#'
#' This dataset includes the topological charge indices descriptors of
#' the 20 amino acids calculated by Dragon (version 5.4)
#' used for scales extraction in this package.
#'
#' @docType data
#' @name AATopoChg
#' @examples
#' data(AATopoChg)
NULL

#' Walk and Path Counts Descriptors for 20 Amino Acids calculated by Dragon
#'
#' This dataset includes the walk and path counts descriptors of
#' the 20 amino acids calculated by Dragon (version 5.4)
#' used for scales extraction in this package.
#'
#' @docType data
#' @name AAWalk
#' @examples
#' data(AAWalk)
NULL

#' WHIM Descriptors for 20 Amino Acids calculated by Dragon
#'
#' This dataset includes the WHIM descriptors of
#' the 20 amino acids calculated by Dragon (version 5.4)
#' used for scales extraction in this package.
#'
#' @docType data
#' @name AAWHIM
#' @examples
#' data(AAWHIM)
NULL



#' BLOSUM45 Matrix for 20 Amino Acids
#'
#' BLOSUM45 Matrix for the 20 amino acids. The matrix was extracted from the
#' \code{Biostrings} package of Bioconductor.
#'
#' @docType data
#' @name AABLOSUM45
#' @examples
#' data(AABLOSUM45)
NULL

#' BLOSUM50 Matrix for 20 Amino Acids
#'
#' BLOSUM50 Matrix for the 20 amino acids. The matrix was extracted from the
#' \code{Biostrings} package of Bioconductor.
#'
#' @docType data
#' @name AABLOSUM50
#' @examples
#' data(AABLOSUM50)
NULL

#' BLOSUM62 Matrix for 20 Amino Acids
#'
#' BLOSUM62 Matrix for the 20 amino acids. The matrix was extracted from the
#' \code{Biostrings} package of Bioconductor.
#'
#' @docType data
#' @name AABLOSUM62
#' @examples
#' data(AABLOSUM62)
NULL

#' BLOSUM80 Matrix for 20 Amino Acids
#'
#' BLOSUM80 Matrix for the 20 amino acids. The matrix was extracted from the
#' \code{Biostrings} package of Bioconductor.
#'
#' @docType data
#' @name AABLOSUM80
#' @examples
#' data(AABLOSUM80)
NULL

#' BLOSUM100 Matrix for 20 Amino Acids
#'
#' BLOSUM100 Matrix for the 20 amino acids. The matrix was extracted from the
#' \code{Biostrings} package of Bioconductor.
#'
#' @docType data
#' @name AABLOSUM100
#' @examples
#' data(AABLOSUM100)
NULL

#' PAM30 Matrix for 20 Amino Acids
#'
#' PAM30 Matrix for the 20 amino acids. The matrix was extracted from the
#' \code{Biostrings} package of Bioconductor.
#'
#' @docType data
#' @name AAPAM30
#' @examples
#' data(AAPAM30)
NULL

#' PAM40 Matrix for 20 Amino Acids
#'
#' PAM40 Matrix for the 20 amino acids. The matrix was extracted from the
#' \code{Biostrings} package of Bioconductor.
#'
#' @docType data
#' @name AAPAM40
#' @examples
#' data(AAPAM40)
NULL

#' PAM70 Matrix for 20 Amino Acids
#'
#' PAM70 Matrix for the 20 amino acids. The matrix was extracted from the
#' \code{Biostrings} package of Bioconductor.
#'
#' @docType data
#' @name AAPAM70
#' @examples
#' data(AAPAM70)
NULL

#' PAM120 Matrix for 20 Amino Acids
#'
#' PAM120 Matrix for the 20 amino acids. The matrix was extracted from the
#' \code{Biostrings} package of Bioconductor.
#'
#' @docType data
#' @name AAPAM120
#' @examples
#' data(AAPAM120)
NULL

#' PAM250 Matrix for 20 Amino Acids
#'
#' PAM250 Matrix for the 20 amino acids. The matrix was extracted from the
#' \code{Biostrings} package of Bioconductor.
#'
#' @docType data
#' @name AAPAM250
#' @examples
#' data(AAPAM250)
NULL
