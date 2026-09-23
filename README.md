# Shiny Aging Murine Multi-Omic Analyzer (Shiny AMMOA) - Local Version


  - [About](#about)
  - [Citation](#citation)
  - [Setup Guide](#setup-guide)
    - [Two Ways to Run Local Shiny AMMOA](#two-ways-to-run-local-shiny-ammoa)
    - [Local Installation Guide](#local-installation-guide)
      - [Install Required Softwares](#install-required-softwares)
      - [Install Required Packages](#install-required-packages)
      - [Download Source Code and Data](#download-source-code-and-data)
      - [Run the App](#run-the-app)
    - [Docker Installation Guide](#docker-installation-guide)
  - [How to Use the App](#how-to-use-the-app)
  - [Troubleshooting](#troubleshooting)


## About 

Shiny Aging Murine Multi-Omic Analyzer (Shiny AMMOA) is a graphical user interface (GUI)-based analytical platform that enables interactive exploration of murine aging-associated bulk transcriptomic, proteomic, and metabolomic datasets.
<br>
<br>

## Citation  

When using this tool, please cite:
> **Shiny AMMOA: an interactive platform for integrative multi-omics analysis of murine aging**  
> Mayuka Ninomiya Kanda  
> bioRxiv 2026.05.18.726091; [DOI link](https://doi.org/10.64898/2026.05.18.726091)

This repository provides the source code and setup guide to start Shiny AMMOA.<br>
<br>


## Setup Guide  

### Two Ways to Run Local Shiny AMMOA

There are two ways to run Shiny AMMOA locally (i.e., in your machine):

1. **Local Installation**  
   Shiny AMMOA runs directly in your local R environment. 
   - *Pros: You can get started relatively quickly without using Docker.*
   - *Cons: This option may be more susceptible to software and package version compatibility issues.*

2. **Docker Installation**  
   Shiny AMMOA runs inside a Docker environment.  
   - *Pros: You can run Shiny AMMOA in the same environment in which the app was implemented and tested, with less concern about package version compatibility.*
   - *Cons: You may need additional computer resources and some technical knowledge to build the Docker image.*

If you are not familiar with Docker, I recommend trying **[Local Installation](#local-installation-guide)** first. If you encounter compatibility or installation issues, you can then try **[Docker Installation](#docker-installation-guide)**.
<br>
<br>

### Local Installation Guide

#### Install Required Softwares
If you do not have the following software installed, please install them. Be sure to install versions compatible with your operating system (Windows / macOS / Linux).

- [R](https://cran.r-project.org/) (ver. ≥ 4.6.0)
- [R Studio](https://posit.co/download/rstudio-desktop) (ver. ≥ 2026.08.2+200)
<br>

#### Install Required Packages
The packages below needs to be installed to run Shiny AMMOA app.

**CRAN Packages**
- shiny (≥ 1.13.0)
- bslib (≥ 0.11.0)
- bsplus (≥ 0.1.5)
- reactable (≥ 0.4.5)
- waiter (≥ 0.2.5.1)
- httr (≥ 1.4.8)
- htmltools (≥ 0.5.9)
- XML (≥ 3.99-0.23)
- ggiraph (≥ 0.9.6)
- tidyverse (≥ 2.0.0)
- RColorBrewer (≥ 1.1-3)
- viridis (≥ 0.6.5)
- png (≥ 0.1-9)
  
**BIoconductor Packages**
- clusterProfiler (≥ 4.20.0)
- MicrobiomeProfiler (≥ 1.18.1)
- enrichplot (≥ 1.32.0)
- org.Mm.eg.db (≥ 3.23.0) 
- pathview (≥ 1.52.0)

The required packages will be installed by running the script below.
To run the code, open RStudio, paste the code into the Console (the panel where you can type commands), and press Enter. Alternatively, you can paste the code into a script file and click “Run”.
The installation process may take anywhere from a few minutes to about 10–15 minutes, depending on your network connection and system performance.
*<u>Note: package installation is only required once; you do not need to run the installation code every time before using Shiny AMMOA.</u>*

```r
# Necessary Packages

# CRAN Packages
cran_pkgs <- c(
  "shiny", "bslib", "bsplus", "reactable", "waiter",
  "httr", "htmltools", "XML",
  "RColorBrewer", "viridis", "ggiraph", "png",
  "tidyverse"
)

# BIoconductor Packages
bioc_pkgs <- c(
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
<br>

#### Download Source Code and Data  
**Step1**.  Go to the Zenodo repository:<br>
[https://doi.org/10.5281/zenodo.22903473](https://doi.org/10.5281/zenodo.22903473)
<br>

**Step 2**. Download the `Shiny_AMMOA_Docker.zip` file and unzip it.
<br>

**Step 3**. Save the `Shiny AMMOA` folder extracted from the ZIP file to any location on your local drive.  
*NOTE: Be sure not to change the folder structure, as this may break the data structure referenced in the app.*
<br>
<br>

#### Run the App  
Open RStudio and run the following commands in the Console:

```r
setwd("your_directory")
shiny::runApp()
```
Replace `your_directory` with the actual path to the extracted folder you saved in **Step 3 above** (e.g., `/Users/your_name/Documents/Shiny_AMMOA`).

<br>

### Docker Installation Guide
If you are new to Docker, please refer to one of the many introductory Docker tutorials available online. (You can also use an AI assistant to help you with basic Docker operations!)

**Step 1**. Download Source Code and Data  
Go to the Zenodo repository:  
[https://doi.org/10.5281/zenodo.22903473](https://doi.org/10.5281/zenodo.22903473)  
And download the `Shiny_AMMOA_Docker.zip` file and unzip it.
<br>

**Step 2**. Rename the Folder  
Rename `Shiny_AMMOA` folder in the unzipped file to `Shiny_AMMOA_local`.  
This folder name is required by the Dockerfile and allows the application files to be correctly referenced when the Docker image is built.
<br>

**Step 3**. Build the Docker Image  
The `Shiny_AMMOA_Docker.zip` file contains a `Dockerfile` and a `docker-compose.yml` file.

For Shiny AMMOA local version 2.0, the Dockerfile uses a [Bioconductor Docker image](#https://hub.docker.com/u/bioconductor) as its base, providing an environment for running Shiny applications with Bioconductor version 3.23. I added commands to the `Dockerfile` to install the required CRAN and Bioconductor packages at the specified versions.

Open a terminal and navigate to the folder containing the `docker-compose.yml` file. Then, run:

```
docker compose up --build
```
Docker Compose will build the Docker image and start the Shiny AMMOA application.

> **Note for Apple Silicon users (M1/M2/M3/M4)**
> If you are using a Mac with an Apple Silicon chip (M1, M2, M3, or M4), you may encounter an issue when building the Docker image with the default Dockerfile.
>Before running `docker compose up --build`, please open the Dockerfile and change the following line:
>
> ```FROM bioconductor/shiny:3.23```
> 
>to:
>
>```FROM --platform=linux/amd64 bioconductor/shiny:3.23```
>
>Save the Dockerfile after making this change, and then proceed with the Docker installation steps above.
>*This change is only required for Apple Silicon Macs.*

<br>

**Step 4**. Run the App  

Once the application has started successfully, open the following address in your web browser:  
http://localhost:3838  
To stop the application, press Ctrl+C in the terminal where Docker Compose is running.
<br>
<br>

## How to Use the App
A user guide for performing bioinformatic analysis is provided within the app. Please refer to the **Userguide** tab.
<br>
<br>

## Troubleshooting  
**Q. My app crashes during KEGG analysis**  
You need an internet connection to run the local version of Shiny AMMOA, as some analyses, especially those involving KEGG data, access online resources during execution.
Some VPNs (Virtual Private Networks) may block access to KEGG from R. If your app crashes during KEGG analysis and your R console shows an error message such as `cannot open the connection to 'https://rest.kegg.jp/link/mmu/pathway'`, you may want to try deactivating your VPN and running the analysis again.  
<br>

**Q. How can I restart the app after a crash?** <br>
If your app crashes due to an unexpected error, including the one described above, you can click the **STOP**🛑 icon in the upper-right corner of the RStudio console panel or press the Esc key to stop the application.  If you are not using R Studio, press Esc key. You can then run shiny::runApp() in the R console to start the app again.
