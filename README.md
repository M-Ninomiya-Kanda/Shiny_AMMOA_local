# Shiny Aging Murine Multi Omic Analyzed (Shiny AMMOA) - Local Version

## About 

Shiny Aging Murine Multi Omic Analyzed (Shiny AMMOA) is a graphical user interface (GUI)-based analytical platform that enables interactive exploration of murine aging-associated bulk transcriptomic, proteomic, and metabolomic datasets.

When using this tool, please cite:
> *Coming Soon*

This repository provides the source code and data required to run the app. The source scripts for data preprocessing and a lightweight web version are available in other repositories on the author’s [GitHub page](https://github.com/M-Ninomiya-Kanda).

A setup guide for running Shiny AMMOA locally is provided below.

## Setting-up Guide

### Install Required Softwares
If you do not have the following software installed, please install them. Be sure to install versions compatible with your operating system (Windows / macOS / Linux).

- [R](https://cran.r-project.org/) (ver. ≥ 4.2.2)
- [R Studio](https://posit.co/download/rstudio-desktop) (ver. ≥ 2026.01.0+392)

### Install Required Packages
The packages below needs to be installed to run Shiny AMMOA app.

**CRAN Packages**
- shiny shiny (≥ 1.13.0)
- bslib (≥ 0.10.0)
- bsplus (≥ 0.1.5)
- reactable (≥ 0.4.5)
- waiter (≥ 0.2.5.1)
- httr (≥ 1.4.8)
- htmltools (≥ 0.5.9)
- data.table (≥ 1.18.2.1)
- XML (≥ 3.99-0.23)
- ggiraph (≥ 0.9.6)
- tidyverse (≥ 2.0.0)
  
**BIoconductor Packages**
- DESeq2 (≥ 1.46.0)
- clusterProfiler (≥ 4.14.0)
- MicrobiomeProfiler (≥ 1.12.0)
- enrichplot (≥ 1.26.0)
- org.Mm.eg.db (≥ 3.20.0) 
- pathview (≥ 1.46.0)

The required packages will be installed by running the script below.
To run the code, open RStudio, paste the code into the Console (the panel where you can type commands), and press Enter. Alternatively, you can paste the code into a script file and click “Run”.
The installation process may take anywhere from a few minutes to about 10–15 minutes, depending on your network connection and system performance.
*<u>Note: package installation is only required once; you do not need to run the installation code every time before using Shiny AMMOA.</u>*

```r
# Necessary Packages

# CRAN Packages
cran_pkgs <- c(
  "shiny", "bslib", "bsplus", "reactable", "waiter",
  "httr", "htmltools", "data.table",
  "XML",
  "RColorBrewer", "viridis", "ggiraph", "png",
  "tidyverse"
)

# BIoconductor Packages
bioc_pkgs <- c(
  "DESeq2",
  "clusterProfiler",
  "MicrobiomeProfiler",
  "enrichplot",
  "org.Mm.eg.db",
  "pathview"
)

# Install missing CRAN packages
cran_missing <- cran_pkgs[!cran_pkgs %in% rownames(installed.packages())]
if (length(cran_missing) > 0) {
  install.packages(cran_missing, repos = "https://cloud.r-project.org")
}

# Install missing Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

bioc_missing <- bioc_pkgs[!bioc_pkgs %in% rownames(installed.packages())]
if (length(bioc_missing) > 0) {
  BiocManager::install(bioc_missing, ask = FALSE, update = FALSE)
}
```

## Download Source Code and Data  
Download this GitHub repository as a ZIP file and unzip it. Save the unzipped folder to any location on your local drive. Be sure not to change the folder structure, as this may break the data structure referenced in the app.

## Run the App
Open RStudio and run the following commands in the Console:

```r
setwd("your_directory")
shiny::runApp()
```
Make sure to replace "your_directory" with the actual path to the extracted folder (e.g., "/Users/your_name/Documents/Shiny_AMMOA_local").

## Use the App
A user guide for performing bioinformatic analysis is provided within the app. Please refer to the **Userguide** tab.
