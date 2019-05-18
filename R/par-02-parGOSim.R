.orgdb <- c(
  "anopheles" = "org.Ag.eg.db",
  "arabidopsis" = "org.At.tair.db",
  "bovine" = "org.Bt.eg.db",
  "canine" = "org.Cf.eg.db",
  "chicken" = "org.Gg.eg.db",
  "chimp" = "org.Pt.eg.db",
  "coelicolor" = "org.Sco.eg.db",
  "ecolik12" = "org.EcK12.eg.db",
  "ecsakai" = "org.EcSakai.eg.db",
  "fly" = "org.Dm.eg.db",
  "human" = "org.Hs.eg.db",
  "malaria" = "org.Pf.plasmo.db",
  "mouse" = "org.Mm.eg.db",
  "pig" = "org.Ss.eg.db",
  "rat" = "org.Rn.eg.db",
  "rhesus" = "org.Mmu.eg.db",
  "worm" = "org.Ce.eg.db",
  "xenopus" = "org.Xl.eg.db",
  "yeast" = "org.Sc.sgd.db",
  "zebrafish" = "org.Dr.eg.db"
)

#' Protein Similarity Calculation based on Gene Ontology (GO) Similarity
#'
#' This function calculates the Gene Ontology (GO) similarity
#' between two groups of GO terms or two Entrez gene IDs.
#'
#' @param id1 Character vector. When length > 1: each element is a GO term;
#' when length = 1: the Entrez Gene ID.
#' @param id2 Character vector. When length > 1: each element is a GO term;
#' when length = 1: the Entrez Gene ID.
#' @param type Input type of id1 and id2, \code{'go'} for GO Terms,
#' \code{"gene"} for gene ID.
#' @param ont Default is \code{"MF"}, can be one of \code{"MF"},
#' \code{"BP"}, or \code{"CC"} subontologies.
#' @param organism Organism name. Default is \code{"human"}, can be one of
#' \code{"anopheles"}, \code{"arabidopsis"}, \code{"bovine"}, \code{"canine"},
#' \code{"chicken"}, \code{"chimp"}, \code{"coelicolor"}, \code{"ecolik12"},
#' \code{"ecsakai"}, \code{"fly"}, \code{"human"}, \code{"malaria"},
#' \code{"mouse"}, \code{"pig"}, \code{"rat"}, \code{"rhesus"},
#' \code{"worm"}, \code{"xenopus"}, \code{"yeast"} or \code{"zebrafish"}.
#' Before specifying the organism, please install the corresponding genome wide
#' \href{https://bioconductor.org/packages/release/BiocViews.html#___OrgDb}{annotation data package}
#' for the selected organism.
#' @param measure Default is \code{"Resnik"}, can be one of
#' \code{"Resnik"}, \code{"Lin"}, \code{"Rel"}, \code{"Jiang"}
#' or \code{"Wang"}.
#' @param combine Default is \code{"BMA"}, can be one of \code{"max"},
#' \code{"average"}, \code{"rcmax"} or \code{"BMA"}
#' for combining semantic similarity scores of multiple GO terms
#' associated with proteins.
#'
#' @return Similarity value.
#'
#' @author Nan Xiao <\url{https://nanx.me}>
#'
#' @seealso See \code{\link{parGOSim}} for protein similarity calculation
#' based on Gene Ontology (GO) semantic similarity.
#' See \code{\link{parSeqSim}} for paralleled protein similarity
#' calculation based on Smith-Waterman local alignment.
#'
#' @export twoGOSim
#'
#' @examples
#' \dontrun{
#'
#' # Be careful when testing this since it involves GO similarity computation
#' # and might produce unpredictable results in some environments
#'
#' library("GOSemSim")
#' library("org.Hs.eg.db")
#'
#' # By GO terms
#' go1 <- c("GO:0004022", "GO:0004024", "GO:0004023")
#' go2 <- c("GO:0009055", "GO:0020037")
#' twoGOSim(go1, go2, type = "go", ont = "MF", measure = "Wang")
#'
#' # By Entrez gene id
#' gene1 <- "1956" # EGFR
#' gene2 <- "2261" # FGFR3
#' twoGOSim(gene1, gene2, type = "gene", ont = "BP", measure = "Lin")
#' }
twoGOSim <- function(
  id1, id2,
  type = c("go", "gene"),
  ont = c("MF", "BP", "CC"), organism = "human",
  measure = "Resnik", combine = "BMA") {
  type <- match.arg(type)
  ont <- match.arg(ont)
  godb <- suppressMessages(GOSemSim::godata(OrgDb = .orgdb[organism], ont = ont, computeIC = TRUE))
  if (type == "go") {
    sim <- GOSemSim::mgoSim(
      id1, id2,
      semData = godb,
      measure = measure, combine = combine
    )
  }
  if (type == "gene") {
    sim <- GOSemSim::geneSim(
      id1, id2,
      semData = godb,
      measure = measure, combine = combine
    )
    if ("geneSim" %in% names(sim)) sim <- sim$geneSim
  }
  sim
}

.goPairSim <- function(twoid, golist, godb, measure, combine) {
  id1 <- twoid[1]
  id2 <- twoid[2]
  if (all(golist[[id1]] == "") | all(golist[[id2]] == "")) {
    sim <- 0L
  } else {
    id1good <- 1:length(golist[[id1]])
    id2good <- 1:length(golist[[id2]])
    gid1 <- as.character(golist[[id1]][id1good])
    gid2 <- as.character(golist[[id2]][id2good])
    res <- try(suppressWarnings(
      GOSemSim::mgoSim(gid1, gid2, semData = godb, measure = measure, combine = combine)
    ),
    silent = TRUE
    )
    if (is.numeric(res)) {
      sim <- res
    } else {
      sim <- 0L
    }
  }
  sim
}

#' Protein Similarity Calculation based on Gene Ontology
#' (GO) Similarity
#'
#' This function calculates protein similarity based on
#' Gene Ontology (GO) similarity.
#'
#' @param golist A list, each component contains
#' a character vector of GO terms or one Entrez Gene ID.
#' @param type Input type for \code{golist}, \code{"go"} for GO Terms,
#' \code{"gene"} for gene ID.
#' @param ont Default is \code{"MF"}, can be one of \code{"MF"},
#' \code{"BP"}, or \code{"CC"} subontologies.
#' @param organism Organism name. Default is \code{"human"}, can be one of
#' \code{"anopheles"}, \code{"arabidopsis"}, \code{"bovine"}, \code{"canine"},
#' \code{"chicken"}, \code{"chimp"}, \code{"coelicolor"}, \code{"ecolik12"},
#' \code{"ecsakai"}, \code{"fly"}, \code{"human"}, \code{"malaria"},
#' \code{"mouse"}, \code{"pig"}, \code{"rat"}, \code{"rhesus"},
#' \code{"worm"}, \code{"xenopus"}, \code{"yeast"} or \code{"zebrafish"}.
#' Before specifying the organism, please install the corresponding genome wide
#' \href{https://bioconductor.org/packages/release/BiocViews.html#___OrgDb}{annotation data package}
#' for the selected organism.
#' @param measure Default is \code{"Resnik"}, can be one of
#' \code{"Resnik"}, \code{"Lin"}, \code{"Rel"}, \code{"Jiang"}
#' or \code{"Wang"}.
#' @param combine Default is \code{"BMA"}, can be one of \code{"max"},
#' \code{"average"}, \code{"rcmax"} or \code{"BMA"}
#' for combining semantic similarity scores of multiple GO terms
#' associated with proteins.
#'
#' @return A \code{n} x \code{n} similarity matrix.
#'
#' @author Nan Xiao <\url{https://nanx.me}>
#'
#' @seealso See \code{\link{twoGOSim}} for calculating the
#' GO semantic similarity between two groups of GO terms or two Entrez gene IDs.
#' See \code{\link{parSeqSim}} for paralleled protein similarity
#' calculation based on Smith-Waterman local alignment.
#'
#' @export parGOSim
#'
#' @examples
#' \dontrun{
#'
#' # Be careful when testing this since it involves GO similarity computation
#' # and might produce unpredictable results in some environments
#'
#' library("GOSemSim")
#' library("org.Hs.eg.db")
#'
#' # By GO Terms
#' # AP4B1
#' go1 <- c(
#'   "GO:0005215", "GO:0005488", "GO:0005515",
#'   "GO:0005625", "GO:0005802", "GO:0005905"
#' )
#' # BCAS2
#' go2 <- c(
#'   "GO:0005515", "GO:0005634", "GO:0005681",
#'   "GO:0008380", "GO:0031202"
#' )
#' # PDE4DIP
#' go3 <- c(
#'   "GO:0003735", "GO:0005622", "GO:0005840",
#'   "GO:0006412"
#' )
#' golist <- list(go1, go2, go3)
#' parGOSim(golist, type = "go", ont = "CC", measure = "Wang")
#'
#' # By Entrez gene id
#' genelist <- list(c("150", "151", "152", "1814", "1815", "1816"))
#' parGOSim(genelist, type = "gene", ont = "BP", measure = "Wang")
#' }
parGOSim <- function(
  golist,
  type = c("go", "gene"),
  ont = c("MF", "BP", "CC"), organism = "human",
  measure = "Resnik", combine = "BMA") {
  type <- match.arg(type)
  ont <- match.arg(ont)
  godb <- suppressMessages(GOSemSim::godata(OrgDb = .orgdb[organism], ont = ont, computeIC = TRUE))
  if (type == "gene") {
    gosimmat <- GOSemSim::mgeneSim(unlist(golist), godb, measure = measure, combine = combine, verbose = FALSE)
  }
  if (type == "go") {
    # generate lower matrix index
    idx <- combn(1:length(golist), 2)
    # input is all pair combination
    gosimlist <- vector("list", ncol(idx))
    for (i in 1:ncol(idx)) {
      gosimlist[[i]] <- .goPairSim(rev(idx[, i]), golist = golist, godb, measure, combine)
    }
    # convert list to matrix
    gosimmat <- matrix(0, length(golist), length(golist))
    for (i in 1:length(gosimlist)) gosimmat[idx[2, i], idx[1, i]] <- gosimlist[[i]]
    gosimmat[upper.tri(gosimmat)] <- t(gosimmat)[upper.tri(t(gosimmat))]
    diag(gosimmat) <- 1
  }
  gosimmat
}
