# ============================================================
# Package Installation Script
# CRAN / Bioconductor
# Package versions can be specified individually
# ============================================================

# ------------------------------------------------------------
# 1. Required package versions
# ------------------------------------------------------------

cran_pkgs <- c(
  shiny        = "1.13.0",
  bslib        = "0.11.0",
  bsplus       = "0.1.5",
  reactable    = "0.4.5",
  waiter       = "0.2.5.1",
  httr         = "1.4.8",
  htmltools    = "0.5.9",
  XML          = "3.99-0.23",
  RColorBrewer = "1.1-3",
  viridis      = "0.6.5",
  ggiraph      = "0.9.6",
  png          = "0.1-8",
  tidyverse    = "2.0.0"
)

bioc_pkgs <- c(
  clusterProfiler   = "4.20.0",
  MicrobiomeProfiler = "1.18.0",
  enrichplot        = "1.32.0",
  org.Mm.eg.db      = "3.23.0",
  pathview          = "1.52.0"
)

# Bioconductor version
bioc_version <- "3.23"


# ------------------------------------------------------------
# 2. Install remotes
# ------------------------------------------------------------

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages(
    "remotes",
    repos = "https://cloud.r-project.org"
  )
}


# ------------------------------------------------------------
# 3. Install CRAN packages with specified versions
# ------------------------------------------------------------

install_cran_versioned <- function(pkgs) {

  for (pkg in names(pkgs)) {

    version <- pkgs[[pkg]]

    installed_version <- if (pkg %in% rownames(installed.packages())) {
      as.character(packageVersion(pkg))
    } else {
      NA_character_
    }

    if (is.na(installed_version)) {

      message(
        "Installing CRAN package: ",
        pkg,
        " (version ", version, ")"
      )

      remotes::install_version(
        package = pkg,
        version = version,
        repos = "https://cloud.r-project.org",
        upgrade = "never"
      )

    } else if (utils::compareVersion(installed_version, version) != 0) {

      message(
        "Updating CRAN package: ",
        pkg,
        " from ",
        installed_version,
        " to ",
        version
      )

      remotes::install_version(
        package = pkg,
        version = version,
        repos = "https://cloud.r-project.org",
        upgrade = "never"
      )

    } else {

      message(
        "Already installed: ",
        pkg,
        " (version ", installed_version, ")"
      )
    }
  }
}


# ------------------------------------------------------------
# 4. Install BiocManager
# ------------------------------------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE)) {

  install.packages(
    "BiocManager",
    repos = "https://cloud.r-project.org"
  )
}


# ------------------------------------------------------------
# 5. Set Bioconductor version
# ------------------------------------------------------------

if (
  !identical(
    as.character(BiocManager::version()),
    bioc_version
  )
) {

  message(
    "Setting Bioconductor version to ",
    bioc_version
  )

  BiocManager::install(
    version = bioc_version,
    ask = FALSE
  )
}


# ------------------------------------------------------------
# 6. Install Bioconductor packages
# ------------------------------------------------------------
#
# BiocManager installs packages according to the selected
# Bioconductor release.
#
# For exact package versions, remotes::install_version()
# is used after selecting the Bioconductor release.
# ------------------------------------------------------------

install_bioc_versioned <- function(pkgs) {

  for (pkg in names(pkgs)) {

    version <- pkgs[[pkg]]

    installed_version <- if (pkg %in% rownames(installed.packages())) {
      as.character(packageVersion(pkg))
    } else {
      NA_character_
    }

    if (is.na(installed_version)) {

      message(
        "Installing Bioconductor package: ",
        pkg,
        " (version ", version, ")"
      )

      BiocManager::install(
        pkg,
        ask = FALSE,
        update = FALSE
      )

    } else if (utils::compareVersion(installed_version, version) != 0) {

      message(
        "Installed version of ",
        pkg,
        " is ",
        installed_version,
        "; requested version is ",
        version
      )

      # Try to install the exact version from Bioconductor archive
      remotes::install_version(
        package = pkg,
        version = version,
        repos = BiocManager::repositories(),
        upgrade = "never"
      )

    } else {

      message(
        "Already installed: ",
        pkg,
        " (version ", installed_version, ")"
      )
    }
  }
}


# ------------------------------------------------------------
# 7. Run installation
# ------------------------------------------------------------

install_cran_versioned(cran_pkgs)

install_bioc_versioned(bioc_pkgs)


# ------------------------------------------------------------
# 8. Verify installed versions
# ------------------------------------------------------------

all_pkgs <- c(cran_pkgs, bioc_pkgs)

installed_versions <- sapply(
  names(all_pkgs),
  function(pkg) {
    if (requireNamespace(pkg, quietly = TRUE)) {
      as.character(packageVersion(pkg))
    } else {
      NA_character_
    }
  }
)

version_check <- data.frame(
  Package = names(all_pkgs),
  Requested = unname(all_pkgs),
  Installed = unname(installed_versions),
  Match = mapply(
    function(requested, installed) {
      !is.na(installed) &&
        utils::compareVersion(installed, requested) == 0
    },
    unname(all_pkgs),
    unname(installed_versions)
  ),
  row.names = NULL
)

print(version_check)
