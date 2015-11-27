# protr

Comprehensive toolkit for generating various numerical representation schemes of protein sequence. The descriptors included in the protr package are extensively utilized in bioinformatics and chemogenomics research.

## Package Description

### Commonly used descriptors

  * Amino acid composition
    * Amino acid composition
    * Dipeptide composition
    * Tripeptide composition

  * Autocorrelation
    * Normalized Moreau-Broto autocorrelation
    * Moran autocorrelation
    * Geary autocorrelation

  * CTD
    * Composition
    * Transition
    * Distribution

  * Conjoint Triad

  * Quasi-sequence-order descriptors
    * Sequence-order-coupling number
    * Quasi-sequence-order descriptors
  
  * Pseudo amino acid composition
    * Pseudo amino acid composition
    * Amphiphilic pseudo amino acid composition

  * Profile-based descriptors
    * Profile-based descriptors derived by PSSM (Position-Specific Scoring Matrix)

### Proteochemometric (PCM) modeling descriptors

  * Scales-based descriptors derived by principal components analysis
    * Scales-based descriptors derived by amino acid properties (AAindex)
    * Scales-based descriptors derived by 20+ classes of 2D and 3D molecular descriptors (Topological, WHIM, VHSE, etc.)
  * Scales-based descriptors derived by factor analysis
  * Scales-based descriptors derived by multidimensional scaling
  * BLOSUM and PAM matrix-derived descriptors

### Similarity Computation

Local and global pairwise sequence alignment for protein sequences:

  * Between two protein sequences
  * Parallelized pairwise similarity calculation with a list of protein sequences

GO semantic similarity measures:

  * Between two groups of GO terms / two Entrez Gene IDs
  * Parallelized pairwise similarity calculation with a list of GO terms / Entrez Gene IDs

### Miscellaneous tools and datasets

  * Retrieve protein sequences from UniProt
  * Read protein sequences in FASTA format
  * Read protein sequences in PDB format
  * Sanity check of the amino acid types appeared in the protein sequences
  * Protein sequence segmentation
  * Auto cross covariance (ACC) for generating scales-based descriptors of the same length
  * 20+ pre-computed 2D and 3D descriptor sets for the 20 amino acids to use with the scales-based descriptors
  * BLOSUM and PAM matrices for the 20 amino acids
  * Meta information of the 20 amino acids

## Web Server

ProtrWeb, the web server built on protr, is located at:

[http://protrweb.scbdd.com](http://protrweb.scbdd.com)

ProtrWeb does not require any knowledge of R programming for the users, it is a user-friendly and one-click-to-go online platform for computing the descriptors presented in the protr package.

## How to cite

Formatted citation:

Xiao, N., Cao, D.-S., Zhu, M.-F., and Xu, Q.-S. (2015). protr/ProtrWeb: R package and web server for generating various numerical representation schemes of protein sequences. _Bioinformatics_ 31, 1857-1859.

BibTeX entry:

```
@article{Xiao2015,
author = {Xiao, N. and Cao, D.-S. and Zhu, M.-F. and Xu, Q.-S.},
doi = {10.1093/bioinformatics/btv042},
issn = {1367-4803},
journal = {Bioinformatics},
number = {11},
pages = {1857--1859},
title = {{protr/ProtrWeb: R package and web server for generating various numerical representation schemes of protein sequences}},
url = {http://bioinformatics.oxfordjournals.org/cgi/doi/10.1093/bioinformatics/btv042},
volume = {31},
year = {2015}
}
```

## Links

  * CRAN: http://cran.r-project.org/web/packages/protr/
  * GitHub: https://github.com/road2stat/protr/
  * Bug report: https://github.com/road2stat/protr/issues/
