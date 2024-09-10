# protr 1.7-4

## Improvements

- Improved dependency error handling for similarity calculation functions
  that use pairwise sequence alignment (#56):
  - Added upfront checks for Biostrings and pwalign dependencies.
  - Now they throw clear errors early if dependency conditions are not met.
  - The new behavior prevents error messages from appearing in results.

# protr 1.7-3

## Improvements

- For sequence similarity calculations using pairwise sequence alignment,
  the Biostrings version is now detected at runtime to determine if the
  pwalign package is needed. This is because the relevant components have been
  moved from Biostrings to pwalign since Bioconductor 3.19 and Biostrings 2.72.0.
  This enhancement ensures that protr works properly regardless of the versions
  of R, Bioconductor, and Biostrings installed (thanks, @ecrespoSSF, #52).

# protr 1.7-2

## Improvements

- Running `citation("protr")` now gives better output with the BibTeX
  citation key. This is improved by adding the `key` argument to the
  `bibentry()` call in `inst/CITATION` (#49).
- Revised `vignette("protr")` to fix typos and grammar issues.
  Updated images to use `knitr::include_graphics()` chunks,
  resolving pkgdown 2.1.0 accessibility hints for missing alt text (#50).

# protr 1.7-1

## New features

- `crossSetSim()` now gains two new arguments `batches` and `verbose`.

  The `batches` argument allows users to split the similarity computations
  into multiple batches, which is useful when dealing with
  a large number of sequences and limited RAM.
  The `verbose` argument enables progress updates during the computation.
  This brings `crossSetSim()` to feature parity with `parSeqSim()`.
  (thanks, @ofleitas, #41)

- A new function `crossSetSimDisk()` has been implemented as a disk-based
  version of `crossSetSim()`.

  This function follows a similar approach to `parSeqSimDisk()`,
  where partial results from each batch are cached on the hard drive and
  merged at the end. This allows for processing larger protein sequence
  sets that may not fit into RAM (#41).

# protr 1.7-0

## New features

- `crossSetSim()` is added for calculating pairwise similarity between two sets
  of protein sequence based on sequence alignment (thanks, @seb-mueller, #34).

# protr 1.6-3

## Bug fixes

- Fixed a minor bug in `extractProtFP()` and `extractProtFPGap()`
  when `index = NULL` (thanks, @fcampelo, #30).

## Improvements

- Added a comment about `system.file()` usage to avoid confusion
  (thanks, @jonalv, #31).
- Replaced previous CI/CD solutions with GitHub Actions workflows.
- Fixed broken or moved links in function documentation and vignettes.
- Replaced the original logo with a new hex sticker logo.

# protr 1.6-2

## Improvements

- Switched to the tidyverse code style.
- Updated GitHub repository links to reflect the handle change.
- Removed AppVeyor CI due to the frequent Bioconductor dependency installation issues.

# protr 1.6-1

## Improvements

- Added a new argument `batches` to `parSeqSim()`.
  The new argument supports breaking down the pairwise similarity computation
  into smaller batches. This is useful when you have a large number of
  protein sequences, enough number of CPU cores, but not enough RAM to
  compute and hold all the pairwise similarities in a single batch.
  Also, use the other new argument `verbose` to track the computation progress.

## New features

- Added a new function `parSeqSimDisk()`.
  Compared to the in-memory version `parSeqSim()`, this new function caches
  the partial results in each batch to the hard drive and merges the results
  together in the end. This could further reduce the memory usage for parallel
  similarity computations involving a large number of protein sequences.

## Bug fixes

- Fixed an issue in `parGOSim()` that will create minor numerical
  inconsistencies in results due to argument matching.

# protr 1.6-0

## Bug fixes

- Updated `twoGOSim()` and `parGOSim()` to use the latest `GOSemSim` API
  for computing GO based semantic similarity.
  Issues in the code examples are also fixed.
  We thank Denisa Duma for the feedback.

# protr 1.5-2

## Bug fixes

- Fixed the API endpoint issue (from HTTP to HTTPS) in `getUniProt()`.

## Improvements

- Added two new parameters `gap.opening` and `gap.extension` to `parSeqSim()`,
  allowing more flexible tuning of the sequence alignment for more types of
  amino acid sequence data. We thank Dr. Maisa Pinheiro for the feedback.
- Added floating TOC and new CSS style in the vignette to improve navigation
  and readability.

# protr 1.5-1

## New features

- Added a new function `removeGaps()` for removing/replacing gaps (`-`) or
  any irregular characters from protein sequences, to make them suitable
  for feature extraction or sequence alignment based similarity computation.
  We thank Dr. Maisa Pinheiro for the feedback.

# protr 1.5-0

## Bug fixes

- Resolved a critical bug due to improper `ifelse` conditioning
  ([3f6e106](https://github.com/nanxstats/protr/commit/3f6e106c93ab9f28c532547f68b3cd9d5cc3d9b4)) for the
  distribution descriptor in CTD. We thank Jielu Yan from
  the University of Macau for kindly reporting this issue.

## Improvements

- General fixes and improvements for the package vignette.

# protr 1.4-2

## Improvements

- Functions are now organized into sections on the documentation website (<https://nanx.me/protr/reference/>).
- Use system font stack instead of Google Fonts in vignettes to avoid pandoc SSL issue.

# protr 1.4-1

## Improvements

- Converted table images to Markdown tables in the vignette.
- Updated the screenshot of protrweb in the vignette.

# protr 1.4-0

## Improvements

- Migrated from Sweave-based PDF vignette to knitr-based HTML vignette.

# protr 1.3-0

## Improvements

- Fix obsolete URLs.
- Better R code formatting.
- Better function documentation and vignette formatting.

# protr 1.2-1

## Improvements

- New documentation website: <https://nanx.me/protr/>.
- Added Windows continuous integration support using AppVeyor.
- Better R file naming scheme.

# protr 1.2-0

## Improvements

- Added continuous integration.
- Code code style improvements.

# protr 1.1-1

## Bug fixes

- Fix URLs that cannot be accessed by `curl -I -L`:
  1. Use <http://protr.org>.
  2. Remove all inaccessible URLs.

# protr 1.1-0

2015-12-28

## Bug fixes

- Bug fix in `extractCTDD()`.

# protr 1.0-1

## Bug fixes

- Improvements for dealing with boundary cases in several functions (thanks for @koefoed's patches).

## Improvements

- Added citation information.

# protr 0.5-1

## Improvements

- Minor improvements and fixes for documentation.

# protr 0.5-0

2014-12-18

## Improvements

- Added functions allowing users to specify their own classification of the amino acid.
- Documentation improvements.
- Other minor improvements.

# protr 0.4-1

## Improvements

- General documentation improvements.

# protr 0.4-0

2014-09-20

## New features

- Added profile-based descriptors derived by PSSM.

# protr 0.3-0

2014-06-20

## Improvements

- Added an example workflow using protr in the vignette.

# protr 0.2-1

## Improvements

- Added a `LICENSE` file according to CRAN policies.

# protr 0.2-0

## New features

- Second release.
- Added Proteochemometric (PCM) Modeling descriptors, parallelized similarity computation derived by protein sequence alignment and Gene Ontology (GO) semantic similarity measures between a list of protein sequences / GO terms / Entrez gene IDs.
- Added miscellaneous tools and datasets.
- Initial version of scales-based descriptors derived by principal components analysis.
- Initial version of scales-based descriptors derived by AA-properties (AAindex).
- Initial version of scales-based descriptors derived by 20+ classes of 2D and 3D molecular descriptors.
- Initial version of scales-based descriptors derived by factor analysis.
- Initial version of scales-based descriptors derived by multidimensional scaling.
- Initial version of BLOSUM and PAM Matrix-Derived Descriptors.
- Initial version of parallelized pairwise similarity calculation with a list of protein sequences.
- Initial version of pairwise semantic similarity calculation with a list of GO terms / Entrez gene IDs.
- Initial version of Auto Cross Covariance (ACC) for generating scales-based descriptors of the same length.
- Introducing ProtWeb, the web service based on protr: <http://protr.org>.

# protr 0.1-0

## New features

- Initial version.
- First version of Amino Acid Composition descriptor.
- First version of Dipeptide Composition descriptor.
- First version of Tripeptide Composition descriptor.
- First version of Normalized Moreau-Broto Autocorrelation descriptor.
- First version of Moran Autocorrelation descriptor.
- First version of Geary Autocorrelation descriptor.
- First version of CTD - Composition descriptor.
- First version of CTD - Transition descriptor.
- First version of CTD - Distribution descriptor.
- First version of Conjoint Triad descriptor.
- First version of Sequence Order Coupling Number descriptor.
- First version of Quasi-Sequence-Order descriptor.
- First version of Pseudo Amino Acid Composition descriptor.
- First version of Amphiphilic Pseudo Amino Acid Composition descriptor.
- First version of `readFASTA()`.
- First version of `getUniProt()`.
- First version of `protcheck()`.
- First version of `protseg()`.
