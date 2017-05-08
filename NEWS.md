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

- Initial version
- First version of Amino Acid Composition descriptor
- First version of Dipeptide Composition descriptor
- First version of Tripeptide Composition descriptor
- First version of Normalized Moreau-Broto Autocorrelation descriptor
- First version of Moran Autocorrelation descriptor
- First version of Geary Autocorrelation descriptor
- First version of CTD - Composition descriptor
- First version of CTD - Transition descriptor
- First version of CTD - Distribution descriptor
- First version of Conjoint Triad descriptor
- First version of Sequence Order Coupling Number descriptor
- First version of Quasi-Sequence-Order descriptor
- First version of Pseudo Amino Acid Composition descriptor
- First version of Amphiphilic Pseudo Amino Acid Composition descriptor
- First version of `readFASTA()`
- First version of `getUniProt()`
- First version of `protcheck()`
- First version of `protseg()`
