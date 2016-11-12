# protr

[![Build Status](https://travis-ci.org/road2stat/protr.svg?branch=master)](https://travis-ci.org/road2stat/protr)
[![CRAN Version](http://www.r-pkg.org/badges/version/protr)](https://cran.r-project.org/package=protr)
[![Downloads from the RStudio CRAN mirror](http://cranlogs.r-pkg.org/badges/protr)](http://cranlogs.r-pkg.org/badges/protr)

Comprehensive toolkit for generating various numerical representation schemes of protein sequence. The descriptors included in the protr package are extensively utilized in bioinformatics and chemogenomics.

## Descriptors List

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

[http://protr.org](http://protr.org)

ProtrWeb does not require any knowledge of R programming for the users, it is a user-friendly and one-click-to-go online platform for computing the descriptors presented in the protr package.

## Paper Citation

Formatted citation:

Nan Xiao, Dong-Sheng Cao, Min-Feng Zhu, and Qing-Song Xu. (2015). protr/ProtrWeb: R package and web server for generating various numerical representation schemes of protein sequences. _Bioinformatics_ 31 (11), 1857-1859.

BibTeX entry:

```
@article{Xiao2015,
author = {Xiao, Nan and Cao, Dong-Sheng and Zhu, Min-Feng and Xu, Qing-Song.},
title = {{protr/ProtrWeb: R package and web server for generating various numerical representation schemes of protein sequences}},
journal = {Bioinformatics},
year = {2015},
volume = {31},
number = {11},
pages = {1857--1859},
doi = {10.1093/bioinformatics/btv042},
issn = {1367-4803},
url = {http://bioinformatics.oxfordjournals.org/content/31/11/1857}
}
```

## Links

  * CRAN: https://cran.r-project.org/package=protr
  * GitHub: https://github.com/road2stat/protr
  * Bug report: https://github.com/road2stat/protr/issues
