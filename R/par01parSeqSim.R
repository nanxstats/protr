.seqPairSim = function (twoid) {
  
  id1 = twoid[1]
  id2 = twoid[2]
  
  if ( protlist[[id1]] == '' | protlist[[id2]] == '' ) {
    
    sim = 0L
    
  } else {
    
    s1  = try(AAString(protlist[[id1]]), silent = TRUE)
    s2  = try(AAString(protlist[[id2]]), silent = TRUE)
    s12 = try(pairwiseAlignment(s1, s2, type = type, substitutionMatrix = submat, scoreOnly = TRUE), silent = TRUE)
    s11 = try(pairwiseAlignment(s1, s1, type = type, substitutionMatrix = submat, scoreOnly = TRUE), silent = TRUE)
    s22 = try(pairwiseAlignment(s2, s2, type = type, substitutionMatrix = submat, scoreOnly = TRUE), silent = TRUE)
    
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
#' @note Note
#' 
#' @keywords Needleman-Wunsch Smith-Waterman local global sequence alignment parallel similarity parSeqSim
#'
#' @aliases parSeqSim
#' 
#' @author Nan Xiao <\url{http://www.road2stat.com}>
#' 
#' @seealso See \code{\link{parGOSim}} for paralleled protein similarity
#' calculation based on Gene Ontology (GO) similarity.
#' 
#' @export parSeqSim
#' 
#' @examples
#' \dontrun{
#' require(Biostrings)
#' require(doMC)  # For Linux/OS X
#' registerDoMC(cores = detectCores())  # For Linux/OS X
#' # require(doParallel)  # For Windows
#' # registerDoParallel(cores = detectCores())  # For Windows
#' parSeqSim(protlist, type = 'local', submat = 'BLOSUM62') }
#' 

parSeqSim = function (protlist, type = 'local', submat = 'BLOSUM62') {
  
  require(Biostrings)
  require(foreach)
  if (.Platform$OS.type == 'windows') {
    require(doParallel)
    registerDoParallel(cores)
  } else {
    require(doMC)
    registerDoMC(cores)
  }
  
  submat = get(submat)
  
  # generate lower matrix index
  idx = combn(1:length(protlist), 2)
  
  # then use foreach parallelization
  # input is all pair combination
  
  seqsimlist = vector('list', ncol(idx))
  
  seqsimlist <- foreach (i = 1:length(seqsimlist), .errorhandling = 'pass') %dopar% {
    xxx <- .seqPairSim(rev(idx[, i]))
  }
  
  # convert list to matrix
  seqsimmat = matrix(0, length(protlist), length(protlist))
  for (i in 1:length(seqsimlist)) seqsimmat[idx[2, i], idx[1, i]] = seqsimlist[i]
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
#' @note Note
#' 
#' @keywords Needleman-Wunsch Smith-Waterman local global sequence alignment parallel similarity parSeqSim
#'
#' @aliases twoSeqSim
#' 
#' @author Nan Xiao <\url{http://www.road2stat.com}>
#' 
#' @seealso See \code{\link{parSeqSim}} for paralleled pairwise 
#' protein similarity calculation based on sequence alignment.
#' 
#' @export twoSeqSim
#' 
#' @examples
#' \dontrun{
#' require(Biostrings)
#' twoSeqSim(seq1, seq2, type = 'local', submat = 'BLOSUM62') }
#' 

twoSeqSim = function (seq1, seq2, type = 'local', submat = 'BLOSUM62') {
  
  submat = get(submat)
  
  # sequence alignment for two protein sequences
  s1  = try(AAString(seq1), silent = TRUE)
  s2  = try(AAString(seq2), silent = TRUE)
  s12 = try(pairwiseAlignment(s1, s2, type = type, substitutionMatrix = submat), silent = TRUE)
  
  return(s12)
  
}
