# protr 1.6-1 (2019-02-24)

## Improvements

- Added a new argument `batches` to `parSeqSim()`. The new argument supports breaking down the pairwise similarity computation into smaller batches. This is useful when you have a large number of protein sequences, enough number of CPU cores, but not enough RAM to compute and hold all the pairwise similarities in a single batch. Also, use the other new argument `verbose` to track the computation progress.

## New Features

- Added a new function `parSeqSimDisk()`. Compared to the in-memory version `parSeqSim()`, this new function caches the partial results in each batch to the hard drive and merges the results together in the end. This could further reduce the memory usage for parallel similarity computations involving a large number of protein sequences.

## Bug Fixes

- Fixed an issue in `parGOSim()` that will create minor numerical inconsistencies in results due to argument matching.

# protr 1.6-0 (2019-02-11)

## Bug Fixes

- Updated `twoGOSim()` and `parGOSim()` to use the latest `GOSemSim` API for computing GO based semantic similarity. Issues in the code examples are also fixed. We thank Denisa Duma for the feedback.

# protr 1.5-2 (2018-11-21)

## Bug Fixes

- Fixed the API endpoint issue (from HTTP to HTTPS) in `getUniProt()`.

## Improvements

- Added two new parameters `gap.opening` and `gap.extension` to `parSeqSim()`, allowing more flexible tuning of the sequence alignment for more types of amino acid sequence data. We thank Dr. Maisa Pinheiro for the feedback.
- Added floating TOC and new CSS style in the vignette to improve navigation and readability.

# protr 1.5-1 (2018-07-12)

## New Features

- Added a new function `removeGaps()` for removing/replacing gaps (`-`) or any irregular characters from protein sequences, to make them suitable for feature extraction or sequence alignment based similarity computation. We thank Dr. Maisa Pinheiro for the feedback.

# protr 1.5-0 (2017-11-17)

## Bug Fixes

- Resolved a critical bug due to improper `ifelse` conditioning ([3f6e106](https://github.com/nanxstats/protr/commit/3f6e106c93ab9f28c532547f68b3cd9d5cc3d9b4)) for the distribution descriptor in CTD. We thank Jielu Yan from the University of Macau for kindly reporting this issue.

## Improvements

- General fixes and improvements for the package vignette.

# protr 1.4-2 (2017-09-28)

## Improvements

- The function list is now organized into sections on the package website (https://nanx.me/protr/reference/).
- Use system font stack instead of Google Fonts in vignettes to avoid pandoc SSL issue.

# protr 1.4-1 (2017-07-08)

## Improvements

- Converted table images to markdown tables in the vignette
- Updated the screenshot of protrweb in the vignette

# protr 1.4-0 (2017-06-06)

## Improvements

- Migrated from Sweave-based PDF vignette to knitr-based HTML vignette

# protr 1.3-0 (2017-05-07)

## Improvements

- Fix obsolete URLs
- Better R code formatting
- Better function documentation and vignette formatting

# protr 1.2-1 (2016-12-29)

## Improvements

- New website: https://nanx.me/protr/
- Added Windows continuous integration support using AppVeyor.
- Better R file naming scheme

# protr 1.2-0 (2016-11-12)

## Improvements

- Added continuous integration
- Code style improvements

# protr 1.1-1 (2015-12-29)

## Bug Fixes

- Fix URLs that cannot be accessed by `curl -I -L`:

    1. Use http://protr.org
    2. Remove all inaccessible URLs

# protr 1.1-0 (2015-12-28)

## Bug Fixes

- Bug fix in `extractCTDD()`

# protr 1.0-1 (2015-11-26)

## Bug Fixes

- Improvements for dealing with boundary cases in several functions (thanks for @koefoed's patches)

## Improvements

- Added citation information

# protr 0.5-1 (2014-12-22)

## Improvements

- Minor improvements and fixes for documentation

# protr 0.5-0 (2014-12-18)

## Improvements

- Added functions allowing users to specify their own classification of the amino acid
- Documentation improvements
- Other minor improvements

# protr 0.4-1 (2014-10-10)

## Improvements

- General documentation improvements

# protr 0.4-0 (2014-09-20)

## New Features

- Added profile-based descriptors derived by PSSM

# protr 0.3-0 (2014-06-20)

## Improvements

- Added example workflow using protr in the vignette

# protr 0.2-1 (2014-01-25)

## Improvements

- Added LICENSE file according to CRAN policies

# protr 0.2-0 (2013-12-10)

## New Features

- second release
- added Proteochemometric (PCM) Modeling descriptors, parallellized similarity computation derived by protein sequence alignment and Gene Ontology (GO) semantic similarity measures between a list of protein sequences / GO terms / Entrez Gene IDs
- added misc tools and datasets
- initial version of Scales-Based Descriptors derived by Principal Components Analysis
- initial version of Scales-Based Descriptors derived by AA-Properties (AAindex)
- initial version of Scales-Based Descriptors derived by 20+ classes of 2D and 3D Molecular Descriptors
- initial version of Scales-Based Descriptors derived by Factor Analysis
- initial version of Scales-Based Descriptors derived by Multidimensional Scaling
- initial version of BLOSUM and PAM Matrix-Derived Descriptors
- initial version of parallelized pairwise similarity calculation with a list of protein sequences
- initial version of pairwise semantic similarity calculation with a list of GO terms / Entrez Gene IDs
- initial version of Auto Cross Covariance (ACC) for generating scales-based descriptors of the same length
- introducing ProtWeb, the web service based on protr: http://protr.org

# protr 0.1-0 (2012-11-18)

## New Features

- initial version
- first version of Amino Acid Composition descriptor
- first version of Dipeptide Composition descriptor
- first version of Tripeptide Composition descriptor
- first version of Normalized Moreau-Broto Autocorrelation descriptor
- first version of Moran Autocorrelation descriptor
- first version of Geary Autocorrelation descriptor
- first version of CTD - Composition descriptor
- first version of CTD - Transition descriptor
- first version of CTD - Distribution descriptor
- first version of Conjoint Triad descriptor
- first version of Sequence Order Coupling Number descriptor
- first version of Quasi-Sequence-Order descriptor
- first version of Pseudo Amino Acid Composition descriptor
- first version of Amphiphilic Pseudo Amino Acid Composition descriptor
- first version of `readFASTA()`
- first version of `getUniProt()`
- first version of `protcheck()`
- first version of `protseg()`
