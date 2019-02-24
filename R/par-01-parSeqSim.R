#' Parallellized Protein Sequence Similarity Calculation based on
#' Sequence Alignment (In-Memory Version)
#'
#' This function implemented the parallellized version for calculating
#' protein sequence similarity based on sequence alignment.
#'
#' @param protlist A length \code{n} list containing \code{n} protein sequences,
#' each component of the list is a character string, storing one protein
#' sequence. Unknown sequences should be represented as \code{""}.
#' @param cores Integer. The number of CPU cores to use for parallel execution,
#' default is \code{2}. Users can use the \code{detectCores()} function
#' in the \code{parallel} package to see how many cores they could use.
#' @param batches Integer. How many batches should we split the pairwise
#' similarity computations into. This is useful when you have a large
#' number of protein sequences, enough number of CPU cores, but not
#' enough RAM to compute and hold all the pairwise similarities
#' in a single batch. Defaults to 1.
#' @param verbose Print the computation progress?
#' @param type Type of alignment, default is \code{"local"},
#' can be \code{"global"} or \code{"local"},
#' where \code{"global"} represents Needleman-Wunsch global alignment;
#' \code{"local"} represents Smith-Waterman local alignment.
#' @param submat Substitution matrix, default is \code{"BLOSUM62"},
#' can be one of \code{"BLOSUM45"}, \code{"BLOSUM50"}, \code{"BLOSUM62"},
#' \code{"BLOSUM80"}, \code{"BLOSUM100"}, \code{"PAM30"},
#' \code{"PAM40"}, \code{"PAM70"}, \code{"PAM120"}, or \code{"PAM250"}.
#' @param gap.opening The cost required to open a gap of any length
#' in the alignment. Defaults to 10.
#' @param gap.extension The cost to extend the length of an existing
#' gap by 1. Defaults to 4.
#'
#' @return A \code{n} x \code{n} similarity matrix.
#'
#' @keywords alignment parallel similarity
#'
#' @aliases parSeqSim
#'
#' @author Nan Xiao <\url{https://nanx.me}>
#'
#' @seealso See \code{\link{parSeqSimDisk}} for the disk-based version.
#'
#' @export parSeqSim
#'
#' @examples
#' \dontrun{
#'
#' # Be careful when testing this since it involves parallelisation
#' # and might produce unpredictable results in some environments
#'
#' library("Biostrings")
#' library("foreach")
#' library("doParallel")
#'
#' s1 <- readFASTA(system.file("protseq/P00750.fasta", package = "protr"))[[1]]
#' s2 <- readFASTA(system.file("protseq/P08218.fasta", package = "protr"))[[1]]
#' s3 <- readFASTA(system.file("protseq/P10323.fasta", package = "protr"))[[1]]
#' s4 <- readFASTA(system.file("protseq/P20160.fasta", package = "protr"))[[1]]
#' s5 <- readFASTA(system.file("protseq/Q9NZP8.fasta", package = "protr"))[[1]]
#' plist <- list(s1, s2, s3, s4, s5)
#' (psimmat <- parSeqSim(plist, cores = 2, type = "local", submat = "BLOSUM62"))
#' }
parSeqSim <- function(
  protlist,
  cores = 2, batches = 1, verbose = FALSE,
  type = "local", submat = "BLOSUM62", gap.opening = 10, gap.extension = 4) {
  doParallel::registerDoParallel(cores)

  # generate lower matrix index
  idx <- combn(1:length(protlist), 2)

  # split index into k batches
  split2 <- function(x, k) split(x, sort(rank(x) %% k))
  idxbatch <- split2(1:ncol(idx), batches)

  # then use foreach parallelization
  # input is all pair combinations (in each batch)
  `%mydopar%` <- foreach::`%dopar%`
  seqsimlist_batch <- vector("list", batches)
  for (k in 1:batches) {
    if (verbose) cat("Starting batch", k, "of", batches, "\n")
    seqsimlist_batch[[k]] <- foreach::foreach(
      i = idxbatch[[k]], .errorhandling = "pass"
    ) %mydopar% {
      tmp <- .seqPairSim(
        rev(idx[, i]), protlist, type, submat, gap.opening, gap.extension
      )
    }
  }

  # merge all batches
  seqsimlist <- as.list(unlist(seqsimlist_batch))

  # convert list to matrix
  seqsimmat <- matrix(0, length(protlist), length(protlist))
  for (i in 1:length(seqsimlist))
    seqsimmat[idx[2, i], idx[1, i]] <- seqsimlist[[i]]
  seqsimmat[upper.tri(seqsimmat)] <- t(seqsimmat)[upper.tri(t(seqsimmat))]
  diag(seqsimmat) <- 1

  seqsimmat
}

#' Parallellized Protein Sequence Similarity Calculation based on
#' Sequence Alignment (Disk-Based Version)
#'
#' This function implemented the parallellized version for calculating
#' protein sequence similarity based on sequence alignment.
#' This version caches the partial results in each batch to the
#' hard drive and merges the results together in the end, which
#' reduces the memory usage.
#'
#' @param protlist A length \code{n} list containing \code{n} protein sequences,
#' each component of the list is a character string, storing one protein
#' sequence. Unknown sequences should be represented as \code{""}.
#' @param cores Integer. The number of CPU cores to use for parallel execution,
#' default is \code{2}. Users can use the \code{detectCores()} function
#' in the \code{parallel} package to see how many cores they could use.
#' @param batches Integer. How many batches should we split the pairwise
#' similarity computations into. This is useful when you have a large
#' number of protein sequences, enough number of CPU cores, but not
#' enough RAM to compute and hold all the pairwise similarities
#' in a single batch. Defaults to 1.
#' @param path Directory for caching the results in each batch.
#' Defaults to the temporary directory.
#' @param verbose Print the computation progress?
#' @param type Type of alignment, default is \code{"local"},
#' can be \code{"global"} or \code{"local"},
#' where \code{"global"} represents Needleman-Wunsch global alignment;
#' \code{"local"} represents Smith-Waterman local alignment.
#' @param submat Substitution matrix, default is \code{"BLOSUM62"},
#' can be one of \code{"BLOSUM45"}, \code{"BLOSUM50"}, \code{"BLOSUM62"},
#' \code{"BLOSUM80"}, \code{"BLOSUM100"}, \code{"PAM30"},
#' \code{"PAM40"}, \code{"PAM70"}, \code{"PAM120"}, or \code{"PAM250"}.
#' @param gap.opening The cost required to open a gap of any length
#' in the alignment. Defaults to 10.
#' @param gap.extension The cost to extend the length of an existing
#' gap by 1. Defaults to 4.
#'
#' @return A \code{n} x \code{n} similarity matrix.
#'
#' @keywords alignment parallel similarity
#'
#' @aliases parSeqSimDisk
#'
#' @seealso See \code{\link{parSeqSim}} for the in-memory version.
#'
#' @author Nan Xiao <\url{https://nanx.me}>
#'
#' @export parSeqSimDisk
#'
#' @examples
#' \dontrun{
#'
#' # Be careful when testing this since it involves parallelisation
#' # and might produce unpredictable results in some environments
#'
#' library("Biostrings")
#' library("foreach")
#' library("doParallel")
#'
#' s1 <- readFASTA(system.file("protseq/P00750.fasta", package = "protr"))[[1]]
#' s2 <- readFASTA(system.file("protseq/P08218.fasta", package = "protr"))[[1]]
#' s3 <- readFASTA(system.file("protseq/P10323.fasta", package = "protr"))[[1]]
#' s4 <- readFASTA(system.file("protseq/P20160.fasta", package = "protr"))[[1]]
#' s5 <- readFASTA(system.file("protseq/Q9NZP8.fasta", package = "protr"))[[1]]
#' set.seed(1010)
#' plist <- as.list(c(s1, s2, s3, s4, s5)[sample(1:5, 100, replace = TRUE)])
#' psimmat <- parSeqSimDisk(
#'   plist,
#'   cores = 2, batches = 10, verbose = TRUE,
#'   type = "local", submat = "BLOSUM62"
#' )
#' }
parSeqSimDisk <- function(
  protlist,
  cores = 2, batches = 1, path = tempdir(), verbose = FALSE,
  type = "local", submat = "BLOSUM62", gap.opening = 10, gap.extension = 4) {
  doParallel::registerDoParallel(cores)

  if (!dir.exists(path)) dir.create(path)

  # generate lower matrix index
  idx <- combn(1:length(protlist), 2)

  # split index into k batches
  split2 <- function(x, k) split(x, sort(rank(x) %% k))
  idxbatch <- split2(1:ncol(idx), batches)

  # then use foreach parallelization
  # input is all pair combinations (in each batch)
  `%mydopar%` <- foreach::`%dopar%`
  for (k in 1:batches) {
    if (verbose) cat("Starting batch", k, "of", batches, "\n")
    seqsimlist_batch_tmp <- foreach::foreach(
      i = idxbatch[[k]], .errorhandling = "pass"
    ) %mydopar% {
      tmp <- .seqPairSim(
        rev(idx[, i]), protlist, type, submat, gap.opening, gap.extension
      )
    }
    # save each batch's results to disk
    saveRDS(seqsimlist_batch_tmp, file = paste0(path, "/protr_batch_", k, ".rds"))
  }

  # read from disk
  seqsimlist_batch <- vector("list", batches)
  for (k in 1:batches) {
    seqsimlist_batch[[k]] <- readRDS(paste0(path, "/protr_batch_", k, ".rds"))
  }

  # merge all batches
  seqsimlist <- as.list(unlist(seqsimlist_batch))

  # convert list to matrix
  seqsimmat <- matrix(0, length(protlist), length(protlist))
  for (i in 1:length(seqsimlist))
    seqsimmat[idx[2, i], idx[1, i]] <- seqsimlist[[i]]
  seqsimmat[upper.tri(seqsimmat)] <- t(seqsimmat)[upper.tri(t(seqsimmat))]
  diag(seqsimmat) <- 1

  seqsimmat
}

#' Protein Sequence Alignment for Two Protein Sequences
#'
#' This function implements the sequence alignment between two protein sequences.
#'
#' @param seq1 Character string, containing one protein sequence.
#' @param seq2 Character string, containing another protein sequence.
#' @param type Type of alignment, default is \code{"local"},
#' could be \code{"global"} or \code{"local"},
#' where \code{"global"} represents Needleman-Wunsch global alignment;
#' \code{"local"} represents Smith-Waterman local alignment.
#' @param submat Substitution matrix, default is \code{"BLOSUM62"},
#' can be one of \code{"BLOSUM45"}, \code{"BLOSUM50"}, \code{"BLOSUM62"},
#' \code{"BLOSUM80"}, \code{"BLOSUM100"}, \code{"PAM30"}, \code{"PAM40"},
#' \code{"PAM70"}, \code{"PAM120"}, or \code{"PAM250"}.
#' @param gap.opening The cost required to open a gap of any length
#' in the alignment. Defaults to 10.
#' @param gap.extension The cost to extend the length of an existing
#' gap by 1. Defaults to 4.
#'
#' @return An \code{Biostrings} object containing the alignment scores
#' and other alignment information.
#'
#' @keywords alignment parallel similarity
#'
#' @aliases twoSeqSim
#'
#' @author Nan Xiao <\url{https://nanx.me}>
#'
#' @seealso See \code{\link{parSeqSim}} for paralleled pairwise
#' protein similarity calculation based on sequence alignment.
#' See \code{\link{twoGOSim}} for calculating the GO semantic
#' similarity between two groups of GO terms or two Entrez gene IDs.
#'
#' @export twoSeqSim
#'
#' @examples
#' \dontrun{
#'
#' # Be careful when testing this since it involves sequence alignment
#' # and might produce unpredictable results in some environments
#' library("Biostrings")
#' s1 <- readFASTA(system.file("protseq/P00750.fasta", package = "protr"))[[1]]
#' s2 <- readFASTA(system.file("protseq/P10323.fasta", package = "protr"))[[1]]
#' seqalign <- twoSeqSim(s1, s2)
#' summary(seqalign)
#' score(seqalign)
#' }
twoSeqSim <- function(
  seq1, seq2, type = "local", submat = "BLOSUM62",
  gap.opening = 10, gap.extension = 4) {

  # sequence alignment for two protein sequences
  s1 <- try(Biostrings::AAString(seq1), silent = TRUE)
  s2 <- try(Biostrings::AAString(seq2), silent = TRUE)
  s12 <- try(Biostrings::pairwiseAlignment(
    s1, s2,
    type = type, substitutionMatrix = submat,
    gapOpening = gap.opening, gapExtension = gap.extension
  ),
  silent = TRUE
  )

  s12
}

.seqPairSim <- function(
  twoid, protlist, type, submat, gap.opening, gap.extension) {
  id1 <- twoid[1]
  id2 <- twoid[2]

  if (protlist[[id1]] == "" |
      protlist[[id2]] == "") {
    sim <- 0L
  } else {
    s1 <- try(Biostrings::AAString(protlist[[id1]]), silent = TRUE)
    s2 <- try(Biostrings::AAString(protlist[[id2]]), silent = TRUE)
    s12 <- try(Biostrings::pairwiseAlignment(
      s1, s2,
      type = type, substitutionMatrix = submat, scoreOnly = TRUE,
      gapOpening = gap.opening, gapExtension = gap.extension
    ),
    silent = TRUE
    )
    s11 <- try(Biostrings::pairwiseAlignment(
      s1, s1,
      type = type, substitutionMatrix = submat, scoreOnly = TRUE,
      gapOpening = gap.opening, gapExtension = gap.extension
    ),
    silent = TRUE
    )
    s22 <- try(Biostrings::pairwiseAlignment(
      s2, s2,
      type = type, substitutionMatrix = submat, scoreOnly = TRUE,
      gapOpening = gap.opening, gapExtension = gap.extension
    ),
    silent = TRUE
    )

    if (is.numeric(s12) == FALSE |
        is.numeric(s11) == FALSE |
        is.numeric(s22) == FALSE) {
      sim <- 0L
    } else if (abs(s11) < .Machine$double.eps |
               abs(s22) < .Machine$double.eps) {
      sim <- 0L
    } else {
      sim <- s12 / sqrt(s11 * s22)
    }
  }

  sim
}
