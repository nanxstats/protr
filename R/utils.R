# These functions:
# - `is_installed()`
# - `get_package_version()`
# - `system_file_cached()`
# were sourced from the shiny package version 1.8.0, available at
# <https://github.com/rstudio/shiny>.
#
# For the original version of these functions, please see:
# <https://github.com/rstudio/shiny/blob/v1.8.0/R/staticimports.R>.
#
# The shiny package is licensed under the GNU General Public License version 3.
# For more details on the license, see
# <https://github.com/rstudio/shiny/blob/main/LICENSE>.

is_installed <- function(pkg, version = NULL) {
  installed <- isNamespaceLoaded(pkg) || nzchar(system_file_cached(package = pkg))

  if (is.null(version)) {
    return(installed)
  }

  if (!is.character(version) && !inherits(version, "numeric_version")) {
    # Avoid https://bugs.r-project.org/show_bug.cgi?id=18548
    alert <- if (identical(Sys.getenv("TESTTHAT"), "true")) stop else warning
    alert("`version` must be a character string or a `package_version` or `numeric_version` object.")

    version <- numeric_version(sprintf("%0.9g", version))
  }

  installed && isTRUE(get_package_version(pkg) >= version)
}

get_package_version <- function(pkg) {
  # `utils::packageVersion()` can be slow, so first try the fast path of
  # checking if the package is already loaded.
  ns <- .getNamespace(pkg)
  if (is.null(ns)) {
    utils::packageVersion(pkg)
  } else {
    as.package_version(ns$.__NAMESPACE__.$spec[["version"]])
  }
}

# A wrapper for `system.file()`, which caches the package path because
# `system.file()` can be slow. If a package is not installed, the result won't
# be cached.
system_file_cached <- local({
  pkg_dir_cache <- character()

  function(..., package = "base") {
    if (!is.null(names(list(...)))) {
      stop("All arguments other than `package` must be unnamed.")
    }

    not_cached <- is.na(match(package, names(pkg_dir_cache)))
    if (not_cached) {
      pkg_dir <- system.file(package = package)
      if (nzchar(pkg_dir)) {
        pkg_dir_cache[[package]] <<- pkg_dir
      }
    } else {
      pkg_dir <- pkg_dir_cache[[package]]
    }

    file.path(pkg_dir, ...)
  }
})

# pairwiseAlignment() has moved from Biostrings to pwalign in Biostrings >= 2.72.0.
# To maximize backward compatibility, we determine which package to use by
# detecting the installation status and of version Biostrings at runtime,
# so that protr can works on any R/Bioconductor/Biostrings version.
is_biostrings_installed <- function() {
  is_installed("Biostrings")
}

is_pwalign_needed <- function() {
  is_installed("Biostrings", version = "2.72.0")
}

is_pwalign_installed <- function() {
  is_installed("pwalign")
}

#' @importFrom utils getFromNamespace
get_pwa <- function(ns) {
  getFromNamespace("pairwiseAlignment", ns = ns)
}

resolve_pwa <- function() {
  if (!is_biostrings_installed()) {
    stop("The package \"Biostrings\" is required. Please install it from Bioconductor.", call. = FALSE)
  }

  if (!is_pwalign_needed()) {
    return(get_pwa("Biostrings"))
  }

  if (!is_pwalign_installed()) {
    stop("The package \"pwalign\" is required. Please install it from Bioconductor.", call. = FALSE)
  }

  get_pwa("pwalign")
}
