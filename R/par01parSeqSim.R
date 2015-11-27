.seqPairSim = function (twoid, protlist = protlist, type = type, submat = submat) {

  id1 = twoid[1]
  id2 = twoid[2]

  if ( protlist[[id1]] == '' | protlist[[id2]] == '' ) {

    sim = 0L

  } else {

    s1  = try(Biostrings::AAString(protlist[[id1]]), silent = TRUE)
    s2  = try(Biostrings::AAString(protlist[[id2]]), silent = TRUE)
    s12 = try(Biostrings::pairwiseAlignment(s1, s2, type = type, substitutionMatrix = submat, scoreOnly = TRUE), silent = TRUE)
    s11 = try(Biostrings::pairwiseAlignment(s1, s1, type = type, substitutionMatrix = submat, scoreOnly = TRUE), silent = TRUE)
    s22 = try(Biostrings::pairwiseAlignment(s2, s2, type = type, substitutionMatrix = submat, scoreOnly = TRUE), silent = TRUE)

    if ( is.numeric(s12) == FALSE | is.numeric(s11) == FALSE | is.numeric(s22) == FALSE ) {
      sim = 0L
    } else if ( abs(s11) < .Machine$double.eps | abs(s22) < .Machine$double.eps ) {
      sim = 0L
    } else {
      sim = s12/sqrt(s11 * s22)
    }

  }

  return(sim)

}

#' Parallellized Protein Sequence Similarity Calculation based on Sequence Alignment
#'
#' Parallellized Protein Sequence Similarity Calculation based on Sequence Alignment
#'
#' This function implemented the parallellized version for calculating
#' protein sequence similarity based on sequence alignment.
#'
#' @param protlist A length \code{n} list containing \code{n} protein sequences,
#' each component of the list is a character string, storing one protein sequence.
#' Unknown sequences should be represented as \code{''}.
#' @param cores Integer. The number of CPU cores to use for parallel execution,
#'        default is \code{2}. Users could use the \code{detectCores()} function
#'        in the \code{parallel} package to see how many cores they could use.
#' @param type Type of alignment, default is \code{'local'},
#' could be \code{'global'} or \code{'local'},
#' where \code{'global'} represents Needleman-Wunsch global alignment;
#' \code{'local'} represents Smith-Waterman local alignment.
#' @param submat Substitution matrix, default is \code{'BLOSUM62'}, could be one of
#' \code{'BLOSUM45'}, \code{'BLOSUM50'}, \code{'BLOSUM62'}, \code{'BLOSUM80'}, \code{'BLOSUM100'},
#' \code{'PAM30'}, \code{'PAM40'}, \code{'PAM70'}, \code{'PAM120'}, \code{'PAM250'}.
#'
#' @return A \code{n} x \code{n} similarity matrix.
#'
#' @keywords alignment parallel similarity
#'
#' @aliases parSeqSim
#'
#' @author Nan Xiao <\url{http://nanx.me}>
#'
#' @seealso See \code{twoSeqSim} for protein sequence alignment
#' for two protein sequences. See \code{\link{parGOSim}} for
#' protein similarity calculation based on
#' Gene Ontology (GO) semantic similarity.
#'
#' @export parSeqSim
#'
#' @examples
#' \donttest{
#' # Be careful when testing this since it involves parallelisation
#' # and might produce unpredictable results in some environments
#'
#' require(Biostrings)
#' require(foreach)
#' require(doParallel)
#'
#' s1 = readFASTA(system.file('protseq/P00750.fasta', package = 'protr'))[[1]]
#' s2 = readFASTA(system.file('protseq/P08218.fasta', package = 'protr'))[[1]]
#' s3 = readFASTA(system.file('protseq/P10323.fasta', package = 'protr'))[[1]]
#' s4 = readFASTA(system.file('protseq/P20160.fasta', package = 'protr'))[[1]]
#' s5 = readFASTA(system.file('protseq/Q9NZP8.fasta', package = 'protr'))[[1]]
#' plist = list(s1, s2, s3, s4, s5)
#' psimmat = parSeqSim(plist, cores = 2, type = 'local', submat = 'BLOSUM62')
#' print(psimmat)}

parSeqSim = function (protlist, cores = 2, type = 'local', submat = 'BLOSUM62') {

  doParallel::registerDoParallel(cores)

  # generate lower matrix index
  idx = combn(1:length(protlist), 2)

  # then use foreach parallelization
  # input is all pair combination

  seqsimlist = vector('list', ncol(idx))

  `%mydopar%` = foreach::`%dopar%`

  seqsimlist <- foreach::foreach(i = 1:length(seqsimlist), .errorhandling = 'pass') %mydopar% {
    tmp <- .seqPairSim(rev(idx[, i]), protlist = protlist, type = type, submat = submat)
  }

  # convert list to matrix
  seqsimmat = matrix(0, length(protlist), length(protlist))
  for (i in 1:length(seqsimlist)) seqsimmat[idx[2, i], idx[1, i]] = seqsimlist[[i]]
  seqsimmat[upper.tri(seqsimmat)] = t(seqsimmat)[upper.tri(t(seqsimmat))]
  diag(seqsimmat) = 1

  return(seqsimmat)

}

#' Protein Sequence Alignment for Two Protein Sequences
#'
#' Protein Sequence Alignment for Two Protein Sequences
#'
#' This function implements the sequence alignment between two protein sequences.
#'
#' @param seq1 A character string, containing one protein sequence.
#' @param seq2 A character string, containing another protein sequence.
#' @param type Type of alignment, default is \code{'local'},
#' could be \code{'global'} or \code{'local'},
#' where \code{'global'} represents Needleman-Wunsch global alignment;
#' \code{'local'} represents Smith-Waterman local alignment.
#' @param submat Substitution matrix, default is \code{'BLOSUM62'}, could be one of
#' \code{'BLOSUM45'}, \code{'BLOSUM50'}, \code{'BLOSUM62'}, \code{'BLOSUM80'}, \code{'BLOSUM100'},
#' \code{'PAM30'}, \code{'PAM40'}, \code{'PAM70'}, \code{'PAM120'}, \code{'PAM250'}.
#'
#' @return An Biostrings object containing the scores and other alignment information.
#'
#' @keywords alignment parallel similarity
#'
#' @aliases twoSeqSim
#'
#' @author Nan Xiao <\url{http://nanx.me}>
#'
#' @seealso See \code{\link{parSeqSim}} for paralleled pairwise
#' protein similarity calculation based on sequence alignment.
#' See \code{\link{twoGOSim}} for calculating the
#' GO semantic similarity between two groups of GO terms or two Entrez gene IDs.
#'
#' @export twoSeqSim
#'
#' @examples
#' \donttest{
#' # Be careful when testing this since it involves sequence alignment
#' # and might produce unpredictable results in some environments
#'
#' require(Biostrings)
#'
#' s1 = readFASTA(system.file('protseq/P00750.fasta', package = 'protr'))[[1]]
#' s2 = readFASTA(system.file('protseq/P10323.fasta', package = 'protr'))[[1]]
#' seqalign = twoSeqSim(s1, s2)
#' summary(seqalign)
#' print(seqalign@@score)}

twoSeqSim = function (seq1, seq2, type = 'local', submat = 'BLOSUM62') {

  # sequence alignment for two protein sequences
  s1  = try(Biostrings::AAString(seq1), silent = TRUE)
  s2  = try(Biostrings::AAString(seq2), silent = TRUE)
  s12 = try(Biostrings::pairwiseAlignment(s1, s2, type = type,
                                          substitutionMatrix = submat),
            silent = TRUE)

  return(s12)

}
