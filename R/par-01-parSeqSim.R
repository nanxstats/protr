#' Parallel Protein Sequence Similarity Calculation Based on
#' Sequence Alignment (In-Memory Version)
#'
#' Parallel calculation of protein sequence similarity based on
#' sequence alignment.
#'
#' @param protlist A length \code{n} list containing \code{n} protein sequences,
#'   each component of the list is a character string, storing one protein
#'   sequence. Unknown sequences should be represented as \code{""}.
#' @param cores Integer. The number of CPU cores to use for parallel execution,
#'   default is \code{2}. Users can use the \code{availableCores()} function
#'   in the parallelly package to see how many cores they could use.
#' @param batches Integer. How many batches should we split the pairwise
#'   similarity computations into. This is useful when you have a large
#'   number of protein sequences, enough number of CPU cores, but not
#'   enough RAM to compute and fit all the pairwise similarities
#'   into a single batch. Defaults to 1.
#' @param verbose Print the computation progress?
#'   Useful when \code{batches > 1}.
#' @param type Type of alignment, default is \code{"local"},
#'   can be \code{"global"} or \code{"local"},
#'   where \code{"global"} represents Needleman-Wunsch global alignment;
#'   \code{"local"} represents Smith-Waterman local alignment.
#' @param submat Substitution matrix, default is \code{"BLOSUM62"},
#'   can be one of \code{"BLOSUM45"}, \code{"BLOSUM50"}, \code{"BLOSUM62"},
#'   \code{"BLOSUM80"}, \code{"BLOSUM100"}, \code{"PAM30"},
#'   \code{"PAM40"}, \code{"PAM70"}, \code{"PAM120"}, or \code{"PAM250"}.
#' @param gap.opening The cost required to open a gap of any length
#'   in the alignment. Defaults to 10.
#' @param gap.extension The cost to extend the length of an existing
#'   gap by 1. Defaults to 4.
#'
#' @return A \code{n} x \code{n} similarity matrix.
#'
#' @author Nan Xiao <\url{https://nanx.me}>
#'
#' @seealso See \code{\link{parSeqSimDisk}} for the disk-based version.
#'
#' @importFrom utils combn
#'
#' @export parSeqSim
#'
#' @examples
#' \dontrun{
#'
#' # Be careful when testing this since it involves parallelization
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
  invisible(resolve_pwa())

  doParallel::registerDoParallel(cores)

  # Generate lower matrix index
  idx <- combn(seq_along(protlist), 2)

  # Split index into k batches
  idxbatch <- split2(seq_len(ncol(idx)), batches)

  # Use foreach parallelization, input is all pair combinations (in each batch).
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

  # Merge all batches
  seqsimlist <- as.list(unlist(seqsimlist_batch))

  # Convert list to matrix
  seqsimmat <- matrix(0, length(protlist), length(protlist))
  for (i in seq_along(seqsimlist)) {
    seqsimmat[idx[2, i], idx[1, i]] <- seqsimlist[[i]]
  }
  seqsimmat[upper.tri(seqsimmat)] <- t(seqsimmat)[upper.tri(t(seqsimmat))]
  diag(seqsimmat) <- 1

  seqsimmat
}

#' Parallel Protein Sequence Similarity Calculation Between Two Sets
#' Based on Sequence Alignment (In-Memory Version)
#'
#' Parallel calculation of protein sequence similarity based on
#' sequence alignment between two sets of protein sequences.
#'
#' @param protlist1 A length \code{n} list containing \code{n} protein sequences,
#'   each component of the list is a character string, storing one protein
#'   sequence. Unknown sequences should be represented as \code{""}.
#' @param protlist2 A length \code{n} list containing \code{m} protein sequences,
#'   each component of the list is a character string, storing one protein
#'   sequence. Unknown sequences should be represented as \code{""}.
#' @param cores Integer. The number of CPU cores to use for parallel execution,
#'   default is \code{2}. Users can use the \code{availableCores()} function
#'   in the parallelly package to see how many cores they could use.
#' @param batches Integer. How many batches should we split the
#'   similarity computations into. This is useful when you have a large
#'   number of protein sequences, enough number of CPU cores, but not
#'   enough RAM to compute and fit all the similarities
#'   into a single batch. Defaults to 1.
#' @param verbose Print the computation progress?
#'   Useful when \code{batches > 1}.
#' @param type Type of alignment, default is \code{"local"},
#'   can be \code{"global"} or \code{"local"},
#'   where \code{"global"} represents Needleman-Wunsch global alignment;
#' \code{"local"} represents Smith-Waterman local alignment.
#' @param submat Substitution matrix, default is \code{"BLOSUM62"},
#'   can be one of \code{"BLOSUM45"}, \code{"BLOSUM50"}, \code{"BLOSUM62"},
#'   \code{"BLOSUM80"}, \code{"BLOSUM100"}, \code{"PAM30"},
#'   \code{"PAM40"}, \code{"PAM70"}, \code{"PAM120"}, or \code{"PAM250"}.
#' @param gap.opening The cost required to open a gap of any length
#'   in the alignment. Defaults to 10.
#' @param gap.extension The cost to extend the length of an existing
#'   gap by 1. Defaults to 4.
#'
#' @return A \code{n} x \code{m} similarity matrix.
#'
#' @author Sebastian Mueller <\url{https://alva-genomics.com}>
#'
#' @importFrom utils combn
#'
#' @export crossSetSim
#'
#' @examples
#' \dontrun{
#'
#' # Be careful when testing this since it involves parallelization
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
#'
#' plist1 <- list(s1 = s1, s2 = s2, s4 = s4)
#' plist2 <- list(s3 = s3, s4_again = s4, s5 = s5, s1_again = s1)
#' psimmat <- crossSetSim(plist1, plist2)
#' colnames(psimmat) <- names(plist1)
#' rownames(psimmat) <- names(plist2)
#' print(psimmat)
#' #                 s1         s2         s4
#' # s3       0.10236985 0.18858241 0.05819984
#' # s4_again 0.04921696 0.12124217 1.00000000
#' # s5       0.03943488 0.06391103 0.05714638
#' # s1_again 1.00000000 0.11825938 0.04921696
#' }
crossSetSim <- function(
    protlist1, protlist2,
    type = "local",
    cores = 2,
    batches = 1,
    verbose = FALSE,
    submat = "BLOSUM62",
    gap.opening = 10,
    gap.extension = 4) {
  invisible(resolve_pwa())

  doParallel::registerDoParallel(cores)

  combinations <- expand.grid(seq_along(protlist1), seq_along(protlist2))

  # Split combinations into k batches
  combinations_batch <- split2(seq_len(nrow(combinations)), batches)

  i <- NULL
  results_batch <- vector("list", batches)
  `%mydopar%` <- foreach::`%dopar%`
  for (k in 1:batches) {
    if (verbose) cat("Starting batch", k, "of", batches, "\n")
    results_batch[[k]] <- foreach::foreach(
      i = combinations_batch[[k]],
      .combine = c,
      .errorhandling = "pass",
      .packages = c("Biostrings")
    ) %mydopar% {
      idx1 <- combinations[i, 1]
      idx2 <- combinations[i, 2]

      .seqPairSim(
        c(idx1, idx2 + length(protlist1)),
        c(protlist1, protlist2),
        type,
        submat,
        gap.opening,
        gap.extension
      )
    }
  }

  # Merge all batches
  results <- as.list(unlist(results_batch))

  matrix(
    results,
    nrow = length(protlist2),
    ncol = length(protlist1),
    byrow = TRUE
  )
}

#' Parallel Protein Sequence Similarity Calculation Based on
#' Sequence Alignment (Disk-Based Version)
#'
#' Parallel calculation of protein sequence similarity based on
#' sequence alignment.
#' This version offloads the partial results in each batch to the
#' hard drive and merges the results together in the end, which
#' reduces the memory usage.
#'
#' @param protlist A length \code{n} list containing \code{n} protein sequences,
#'   each component of the list is a character string, storing one protein
#'   sequence. Unknown sequences should be represented as \code{""}.
#' @param cores Integer. The number of CPU cores to use for parallel execution,
#'   default is \code{2}. Users can use the \code{availableCores()} function
#'   in the parallelly package to see how many cores they could use.
#' @param batches Integer. How many batches should we split the pairwise
#'   similarity computations into. This is useful when you have a large
#'   number of protein sequences, enough number of CPU cores, but not
#'   enough RAM to compute and fit all the pairwise similarities
#'   into a single batch. Defaults to 1.
#' @param path Directory for caching the results in each batch.
#'   Defaults to the temporary directory.
#' @param verbose Print the computation progress?
#'   Useful when \code{batches > 1}.
#' @param type Type of alignment, default is \code{"local"},
#'   can be \code{"global"} or \code{"local"},
#'   where \code{"global"} represents Needleman-Wunsch global alignment;
#'   \code{"local"} represents Smith-Waterman local alignment.
#' @param submat Substitution matrix, default is \code{"BLOSUM62"},
#'   can be one of \code{"BLOSUM45"}, \code{"BLOSUM50"}, \code{"BLOSUM62"},
#'   \code{"BLOSUM80"}, \code{"BLOSUM100"}, \code{"PAM30"},
#'   \code{"PAM40"}, \code{"PAM70"}, \code{"PAM120"}, or \code{"PAM250"}.
#' @param gap.opening The cost required to open a gap of any length
#'   in the alignment. Defaults to 10.
#' @param gap.extension The cost to extend the length of an existing
#'   gap by 1. Defaults to 4.
#'
#' @return A \code{n} x \code{n} similarity matrix.
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
#' # Be careful when testing this since it involves parallelization
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
  invisible(resolve_pwa())

  doParallel::registerDoParallel(cores)

  if (!dir.exists(path)) dir.create(path)

  # Generate lower matrix index
  idx <- combn(seq_along(protlist), 2)

  # Split index into k batches
  idxbatch <- split2(seq_len(ncol(idx)), batches)

  # Use foreach parallelization, input is all pair combinations (in each batch).
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
    # Save each batch's results to disk
    saveRDS(seqsimlist_batch_tmp, file = paste0(path, "/protr_batch_", k, ".rds"))
  }

  # Read from disk
  seqsimlist_batch <- vector("list", batches)
  for (k in 1:batches) {
    seqsimlist_batch[[k]] <- readRDS(paste0(path, "/protr_batch_", k, ".rds"))
  }

  # Merge all batches
  seqsimlist <- as.list(unlist(seqsimlist_batch))

  # Convert list to matrix
  seqsimmat <- matrix(0, length(protlist), length(protlist))
  for (i in seq_along(seqsimlist)) {
    seqsimmat[idx[2, i], idx[1, i]] <- seqsimlist[[i]]
  }
  seqsimmat[upper.tri(seqsimmat)] <- t(seqsimmat)[upper.tri(t(seqsimmat))]
  diag(seqsimmat) <- 1

  seqsimmat
}

#' Parallel Protein Sequence Similarity Calculation Between Two Sets
#' Based on Sequence Alignment (Disk-Based Version)
#'
#' Parallel calculation of protein sequence similarity based on
#' sequence alignment between two sets of protein sequences.
#' This version offloads the partial results in each batch to the
#' hard drive and merges the results together in the end, which
#' reduces the memory usage.
#'
#' @param protlist1 A length \code{n} list containing \code{n} protein sequences,
#'   each component of the list is a character string, storing one protein
#'   sequence. Unknown sequences should be represented as \code{""}.
#' @param protlist2 A length \code{n} list containing \code{m} protein sequences,
#'   each component of the list is a character string, storing one protein
#'   sequence. Unknown sequences should be represented as \code{""}.
#' @param cores Integer. The number of CPU cores to use for parallel execution,
#'   default is \code{2}. Users can use the \code{availableCores()} function
#'   in the parallelly package to see how many cores they could use.
#' @param batches Integer. How many batches should we split the pairwise
#'   similarity computations into. This is useful when you have a large
#'   number of protein sequences, enough number of CPU cores, but not
#'   enough RAM to compute and fit all the pairwise similarities
#'   into a single batch. Defaults to 1.
#' @param path Directory for caching the results in each batch.
#'   Defaults to the temporary directory.
#' @param verbose Print the computation progress?
#'   Useful when \code{batches > 1}.
#' @param type Type of alignment, default is \code{"local"},
#'   can be \code{"global"} or \code{"local"},
#'   where \code{"global"} represents Needleman-Wunsch global alignment;
#'   \code{"local"} represents Smith-Waterman local alignment.
#' @param submat Substitution matrix, default is \code{"BLOSUM62"},
#'   can be one of \code{"BLOSUM45"}, \code{"BLOSUM50"}, \code{"BLOSUM62"},
#'   \code{"BLOSUM80"}, \code{"BLOSUM100"}, \code{"PAM30"},
#'   \code{"PAM40"}, \code{"PAM70"}, \code{"PAM120"}, or \code{"PAM250"}.
#' @param gap.opening The cost required to open a gap of any length
#'   in the alignment. Defaults to 10.
#' @param gap.extension The cost to extend the length of an existing
#'   gap by 1. Defaults to 4.
#'
#' @return A \code{n} x \code{m} similarity matrix.
#'
#' @seealso See \code{\link{crossSetSim}} for the in-memory version.
#'
#' @author Nan Xiao <\url{https://nanx.me}>
#'
#' @export crossSetSimDisk
#'
#' @examples
#' \dontrun{
#'
#' # Be careful when testing this since it involves parallelization
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
#'
#' set.seed(1010)
#' plist1 <- as.list(c(s1, s2, s3, s4, s5)[sample(1:5, 100, replace = TRUE)])
#' plist2 <- as.list(c(s1, s2, s3, s4, s5)[sample(1:5, 100, replace = TRUE)])
#' psimmat <- crossSetSimDisk(
#'   plist1, plist2,
#'   cores = 2, batches = 10, verbose = TRUE,
#'   type = "local", submat = "BLOSUM62"
#' )
#' }
crossSetSimDisk <- function(
    protlist1, protlist2,
    cores = 2, batches = 1, path = tempdir(), verbose = FALSE,
    type = "local", submat = "BLOSUM62", gap.opening = 10, gap.extension = 4) {
  invisible(resolve_pwa())

  doParallel::registerDoParallel(cores)

  if (!dir.exists(path)) dir.create(path)

  combinations <- expand.grid(seq_along(protlist1), seq_along(protlist2))

  # Split combinations into batches
  combinations_batch <- split2(seq_len(nrow(combinations)), batches)

  i <- NULL
  `%mydopar%` <- foreach::`%dopar%`
  for (k in 1:batches) {
    if (verbose) cat("Starting batch", k, "of", batches, "\n")
    results_batch_tmp <- foreach::foreach(
      i = combinations_batch[[k]],
      .errorhandling = "pass"
    ) %mydopar% {
      idx1 <- combinations[i, 1]
      idx2 <- combinations[i, 2]

      .seqPairSim(
        c(idx1, idx2 + length(protlist1)),
        c(protlist1, protlist2),
        type,
        submat,
        gap.opening,
        gap.extension
      )
    }
    # Save each batch's results to disk
    saveRDS(results_batch_tmp, file = paste0(path, "/crossSetSim_batch_", k, ".rds"))
  }

  # Read from disk
  results_batch <- vector("list", batches)
  for (k in 1:batches) {
    results_batch[[k]] <- readRDS(paste0(path, "/crossSetSim_batch_", k, ".rds"))
  }

  # Merge all batches
  results <- as.list(unlist(results_batch))

  matrix(
    results,
    nrow = length(protlist2),
    ncol = length(protlist1),
    byrow = TRUE
  )
}

#' Protein Sequence Alignment for Two Protein Sequences
#'
#' Sequence alignment between two protein sequences.
#'
#' @param seq1 Character string, containing one protein sequence.
#' @param seq2 Character string, containing another protein sequence.
#' @param type Type of alignment, default is \code{"local"},
#'   could be \code{"global"} or \code{"local"},
#'   where \code{"global"} represents Needleman-Wunsch global alignment;
#'   \code{"local"} represents Smith-Waterman local alignment.
#' @param submat Substitution matrix, default is \code{"BLOSUM62"},
#'   can be one of \code{"BLOSUM45"}, \code{"BLOSUM50"}, \code{"BLOSUM62"},
#'   \code{"BLOSUM80"}, \code{"BLOSUM100"}, \code{"PAM30"}, \code{"PAM40"},
#'   \code{"PAM70"}, \code{"PAM120"}, or \code{"PAM250"}.
#' @param gap.opening The cost required to open a gap of any length
#'   in the alignment. Defaults to 10.
#' @param gap.extension The cost to extend the length of an existing
#'   gap by 1. Defaults to 4.
#'
#' @return A \code{Biostrings} object containing the alignment scores
#'   and other alignment information.
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
#'
#' library("Biostrings")
#'
#' s1 <- readFASTA(system.file("protseq/P00750.fasta", package = "protr"))[[1]]
#' s2 <- readFASTA(system.file("protseq/P10323.fasta", package = "protr"))[[1]]
#' seqalign <- twoSeqSim(s1, s2)
#' summary(seqalign)
#' score(seqalign)
#' }
twoSeqSim <- function(
    seq1, seq2, type = "local", submat = "BLOSUM62",
    gap.opening = 10, gap.extension = 4) {
  pwa <- resolve_pwa()

  # Sequence alignment for two protein sequences
  s1 <- try(Biostrings::AAString(seq1), silent = TRUE)
  s2 <- try(Biostrings::AAString(seq2), silent = TRUE)
  s12 <- try(
    pwa(
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
  pwa <- resolve_pwa()

  id1 <- twoid[1]
  id2 <- twoid[2]

  if (protlist[[id1]] == "" |
    protlist[[id2]] == "") {
    sim <- 0L
  } else {
    s1 <- try(Biostrings::AAString(protlist[[id1]]), silent = TRUE)
    s2 <- try(Biostrings::AAString(protlist[[id2]]), silent = TRUE)
    s12 <- try(
      pwa(
        s1, s2,
        type = type, substitutionMatrix = submat, scoreOnly = TRUE,
        gapOpening = gap.opening, gapExtension = gap.extension
      ),
      silent = TRUE
    )
    s11 <- try(
      pwa(
        s1, s1,
        type = type, substitutionMatrix = submat, scoreOnly = TRUE,
        gapOpening = gap.opening, gapExtension = gap.extension
      ),
      silent = TRUE
    )
    s22 <- try(
      pwa(
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

split2 <- function(x, k) split(x, sort(rank(x) %% k))
