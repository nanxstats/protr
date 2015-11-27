#' Compute PSSM (Position-Specific Scoring Matrix) for given protein sequence
#'
#' Compute PSSM (Position-Specific Scoring Matrix) for given protein sequence
#'
#' This function calculates the PSSM (Position-Specific Scoring Matrix) derived
#' by PSI-Blast for given protein sequence or peptides.
#' For given protein sequences or peptides, PSSM represents the
#' log-likelihood of the substitution of the 20 types of amino acids at that
#' position in the sequence. Note that the output value is not normalized.
#'
#' @param seq Character vector, as the input protein sequence.
#' @param start.pos Optional integer denoting the start position of the
#' fragment window. Default is \code{1},
#' i.e. the first amino acid of the given sequence.
#' @param end.pos Optional integer denoting the end position of the
#' fragment window. Default is \code{nchar(seq)},
#' i.e. the last amino acid of the given sequence.
#' @param psiblast.path Character string indicating the path of the
#' \code{psiblast} program.
#' If NCBI Blast+ was previously installed in the operation system,
#' the path will be automatically detected.
#' @param makeblastdb.path Character string indicating the path of the
#' \code{makeblastdb} program.
#' If NCBI Blast+ was previously installed in the system,
#' the path will be automatically detected.
#' @param database.path Character string indicating the path of
#' a reference database (a FASTA file).
#' @param iter Number of iterations to perform for PSI-Blast.
#' @param silent Logical. Whether the PSI-Blast running output
#' should be shown or not (May not work on some Windows versions and
#' PSI-Blast versions), default is \code{TRUE}.
#'
#' @param evalue Expectation value (E) threshold for saving hits.
#' Default is \code{10}.
#' @param word.size Word size for wordfinder algorithm. An integer >= 2.
#' @param gapopen Integer. Cost to open a gap.
#' @param gapextend Integer. Cost to extend a gap.
#' @param matrix Character string. The scoring matrix name
#' (default is \code{'BLOSUM62'}).
#' @param threshold Minimum word score such that the word is added to
#' the BLAST lookup table. A real value >= 0.
#' @param seg Character string. Filter query sequence with SEG (\code{'yes'},
#' \code{'window locut hicut'}, or \code{'no'} to disable) Default is \code{'no'}.
#' @param soft.masking Logical. Apply filtering locations as soft masks?
#' Default is \code{FALSE}.
#' @param culling.limit An integer >= 0. If the query range of a hit is
#' enveloped by that of at least this many higher-scoring hits,
#' delete the hit. Incompatible with \code{best.hit.overhang} and
#' \code{best_hit_score_edge}.
#' @param best.hit.overhang Best Hit algorithm overhang value
#' (A real value >= 0 and =< 0.5, recommended value: 0.1).
#' Incompatible with \code{culling_limit}.
#' @param best.hit.score.edge Best Hit algorithm score edge value
#' (A real value >=0 and =< 0.5, recommended value: 0.1).
#' Incompatible with \code{culling_limit}.
#' @param xdrop.ungap X-dropoff value (in bits) for ungapped extensions.
#' @param xdrop.gap X-dropoff value (in bits) for preliminary gapped extensions.
#' @param xdrop.gap.final X-dropoff value (in bits) for final gapped alignment.
#' @param window.size An integer >= 0. Multiple hits window size,
#' To specify 1-hit algorithm, use \code{0}.
#' @param gap.trigger Number of bits to trigger gapping. Default is \code{22}.
#' @param num.threads Integer. Number of threads (CPUs) to use in the
#' BLAST search. Default is \code{1}.
#' @param pseudocount Integer. Pseudo-count value used when constructing PSSM.
#' Default is \code{0}.
#' @param inclusion.ethresh E-value inclusion threshold for pairwise alignments.
#' Default is \code{0.002}.
#'
#' @return The original PSSM, a numeric matrix which has
#' \code{end.pos - start.pos + 1} columns and \code{20} named rows.
#'
#' @note
#' The function requires the \code{makeblastdb} and \code{psiblast} programs
#' to be properly installed in the operation system or
#' their paths provided.
#'
#' The two command-line programs are included in the NCBI-BLAST+
#' software package. To install NCBI Blast+, just open the NCBI FTP site
#' using web browser or FTP software:
#' \url{ftp://anonymous@@ftp.ncbi.nlm.nih.gov:21/blast/executables/blast+/LATEST/}
#' then download the executable version of BLAST+ according to your
#' operation system, and compile or install the downloaded
#' source code or executable program.
#'
#' Ubuntu/Debian users can directly use the command
#' \code{sudo apt-get install ncbi-blast+} to install NCBI Blast+.
#' For OS X users, download \code{ncbi-blast- ... .dmg} then install.
#' For Windows users, download \code{ncbi-blast- ... .exe} then install.
#'
#' @seealso \link{extractPSSMFeature} \link{extractPSSMAcc}
#'
#' @keywords extract PSSM Blast Alignment
#'
#' @aliases extractPSSM
#'
#' @author Nan Xiao <\url{http://nanx.me}>
#'
#' @export extractPSSM
#'
#' @references
#' Altschul, Stephen F., et al.
#' "Gapped BLAST and PSI-BLAST: a new generation of protein database search programs."
#' \emph{Nucleic acids research} 25.17 (1997): 3389--3402.
#'
#' Ye, Xugang, Guoli Wang, and Stephen F. Altschul.
#' "An assessment of substitution scores for protein profile-profile comparison."
#' \emph{Bioinformatics} 27.24 (2011): 3356--3363.
#'
#' Rangwala, Huzefa, and George Karypis.
#' "Profile-based direct kernels for remote homology detection and fold recognition."
#' \emph{Bioinformatics} 21.23 (2005): 4239--4247.
#'
#' @examples
#' if (Sys.which('makeblastdb') == '' | Sys.which('psiblast') == '') {
#'   cat('Could not find makeblastdb or psiblast. Please install NCBI Blast+ first.')
#' } else {
#'   x = readFASTA(system.file('protseq/P00750.fasta', package = 'protr'))[[1]]
#'   dbpath = tempfile('tempdb', fileext = '.fasta')
#'   invisible(file.copy(from = system.file('protseq/Plasminogen.fasta',
#'                                          package = 'protr'), to = dbpath))
#'   pssmmat = extractPSSM(seq = x, database.path = dbpath)
#'   dim(pssmmat)  # 20 x 562 (P00750: length 562, 20 Amino Acids)
#' }

extractPSSM = function(seq, start.pos = 1L, end.pos = nchar(seq),
                       psiblast.path = NULL, makeblastdb.path = NULL,
                       database.path = NULL, iter = 5, silent = TRUE,
                       evalue = 10L, word.size = NULL,
                       gapopen = NULL, gapextend = NULL,
                       matrix = 'BLOSUM62', threshold = NULL,
                       seg = 'no', soft.masking = FALSE,
                       culling.limit = NULL, best.hit.overhang = NULL,
                       best.hit.score.edge = NULL,
                       xdrop.ungap = NULL, xdrop.gap = NULL,
                       xdrop.gap.final = NULL,
                       window.size = NULL, gap.trigger = 22L,
                       num.threads = 1L, pseudocount = 0L,
                       inclusion.ethresh = 0.002) {

  if (Sys.which('makeblastdb') == '' & is.null(makeblastdb.path))
    stop('Please install makeblastdb (included in the NCBI BLAST+) or specify makeblastdb.path.')

  if (Sys.which('psiblast') == '' & is.null(psiblast.path))
    stop('Please install psiblast (included in the NCBI BLAST+) or specify psiblast.path.')

  makeblastdb.path = if (!is.null(makeblastdb.path)) makeblastdb.path else Sys.which('makeblastdb')
  psiblast.path = if (!is.null(psiblast.path)) psiblast.path else Sys.which('psiblast')

  if (is.null(database.path)) stop('Must specify the database (a FASTA file) path')

  N = end.pos - start.pos + 1L

  # Prepare data for PSI-Blast

  tmpdb = tempfile('protrPSIBlastDB')

  cmddb = paste0(shQuote(makeblastdb.path), ' -dbtype prot -in ',
                 shQuote(database.path), ' -out ', shQuote(tmpdb))

  if (silent == TRUE) system(cmddb, ignore.stdout = TRUE) else system(cmddb)

  # Basic parameters for PSI-Blast

  tmp = tempfile('protrPSIBlast')
  queryFasta = paste0(tmp, '.fasta')
  PSSMfile = paste0(tmp, '.pssm')
  querySeq = Biostrings::AAStringSet(seq)
  Biostrings::writeXStringSet(querySeq, queryFasta)

  # Additional parameters for PSI-Blast

  if (!is.null(evalue)) {
    if (evalue <= 0L) {
      stop('evalue must be > 0')
    }
  }

  evalue.arg = ifelse(evalue == 10L, ' ', paste0(' -evalue ', evalue))

  if (!is.null(word.size)) {
    if (word.size < 2L | word.size >= 6L) {
      stop('word.size must be an integer >= 2 and < 6')
    }
  }

  word.size.arg = ifelse(is.null(word.size), ' ', paste0(' -word_size ', word.size))

  gapopen.arg = ifelse(is.null(gapopen), ' ', paste0(' -gapopen ', gapopen))
  gapextend.arg = ifelse(is.null(gapextend), ' ', paste0(' -gapextend ', gapextend))

  matrix.arg = ifelse(matrix == 'BLOSUM62', ' ', paste0(' -matrix ', matrix))

  if (!is.null(threshold)) {
    if (threshold < 0L) {
      stop('threshold must be a real number >= 0')
    }
  }

  threshold.arg = ifelse(is.null(threshold), ' ', paste0(' -threshold ', threshold))

  seg.arg = ifelse(seg == 'no', ' ', paste0(' -seg ', seg))
  soft.masking.arg = ifelse(soft.masking == FALSE, ' ', ' -soft_masking true')

  if (!is.null(culling.limit)) {
    if (culling.limit < 0L) {
      stop('culling.limit must be an integer >= 0')
    }
  }

  culling.limit.arg = ifelse(is.null(culling.limit), ' ', paste0(' -culling_limit ', culling.limit))

  if (!is.null(best.hit.overhang)) {
    if (best.hit.overhang < 0L | best.hit.overhang > 1.5) {
      stop('best.hit.overhang must be a real number >= 0 and =< 0.5')
    }
  }

  best.hit.overhang.arg = ifelse(is.null(best.hit.overhang), ' ', paste0(' -best_hit_overhang ', best.hit.overhang))

  if (!is.null(best.hit.score.edge)) {
    if (best.hit.score.edge < 0L | best.hit.score.edge > 1.5) {
      stop('best.hit.score.edge must be a real number >= 0 and =< 0.5')
    }
  }

  best.hit.score.edge.arg = ifelse(is.null(best.hit.score.edge), ' ', paste0(' -best_hit_score_edge ', best.hit.score.edge))

  xdrop.ungap.arg = ifelse(is.null(xdrop.ungap), ' ', paste0(' -xdrop_ungap ', xdrop.ungap))
  xdrop.gap.arg = ifelse(is.null(xdrop.gap), ' ', paste0(' -xdrop_gap ', xdrop.gap))
  xdrop.gap.final.arg = ifelse(is.null(xdrop.gap.final), ' ', paste0(' -xdrop_gap_final ', xdrop.gap.final))
  window.size.arg = ifelse(is.null(window.size), ' ', paste0(' -window_size ', window.size))
  gap.trigger.arg = ifelse(gap.trigger == 22L, ' ', paste0(' -gap_trigger ', gap.trigger))
  num.threads.arg = ifelse(num.threads == 1L, ' ', paste0(' -num_threads ', num.threads))
  pseudocount.arg = ifelse(pseudocount == 0L, ' ', paste0(' -pseudocount ', pseudocount))
  inclusion.ethresh.arg = ifelse(inclusion.ethresh == 0.002, ' ', paste0(' -inclusion_ethresh ', inclusion.ethresh))

  # Run PSI-Blast

  cmdpsi = paste(
    paste0(shQuote(psiblast.path),
           ' -comp_based_stats 1 -db ', shQuote(tmpdb),
           ' -query ', shQuote(queryFasta), ' -num_iterations ', iter,
           ' -out_ascii_pssm ', shQuote(PSSMfile)),
    paste0(evalue.arg, word.size.arg, gapopen.arg, gapextend.arg,
           matrix.arg, threshold.arg, seg.arg, soft.masking.arg,
           culling.limit.arg, best.hit.overhang.arg, best.hit.score.edge.arg,
           xdrop.ungap.arg, xdrop.gap.arg, xdrop.gap.final.arg, window.size.arg,
           gap.trigger.arg, num.threads.arg, pseudocount.arg,
           inclusion.ethresh.arg))

  if (silent == TRUE) system(cmdpsi, ignore.stdout = TRUE) else system(cmdpsi)

  PSSMraw = readLines(PSSMfile)

  PSSMClean = function(x) {
    y = unlist(strsplit(PSSMraw[x + 3], split = '\\s+'))
    return(-as.numeric(y[4L:23L]))
  }

  PSSMmat = sapply(start.pos:end.pos, PSSMClean)

  AADict = c('A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
             'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V')
  rownames(PSSMmat) = AADict

  return(PSSMmat)

}
