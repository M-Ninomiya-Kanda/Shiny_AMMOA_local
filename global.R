library(shiny)
library(bslib)
library(bsplus)
library(reactable)
library(waiter)

library(httr)
library(htmltools)
library(XML)

library(RColorBrewer)
library(viridis)
library(ggiraph)
library(png)

library(tidyverse)

# load supporting file (mostly in-house functions)
source("shiny_multiomic_source.R")

# Load Necessary Data -----

rna_metadata <- readRDS("data/bulkTMS_metadata_exported_2025_12_17.rds")
gene_ids <- readRDS("data/bulkTMS_ENTREZID_list_standalone_2025_12_28.rds")

# result of proteomic estimation result from each tissue (CellRep 2023)
dat_prot <- readRDS("data/CellRep2023_tissue_res_2025-11-07.rds")

# result of metabolomic estimation result from each tissue (CellMet 2025)
dat_metab <- readRDS("data/CellMet2025_metabolites_adjustedLFC_2025-11-25.rds")
# compound IDs annotations
dat_compound_id <- read.csv("data/compound_mapping_res_2025-11-25.csv")

# load KEGG pathway IDs
kegg_ids <- readRDS("data/KEGG_pathway_ID_list_2025_12_17_exported.rds")


# Define Objects Required for Dynamic UI Change -----
# change tissue list in response to datasource (RNA/Protein)
choices_tissue <- list(
  rna = c("BAT", "Bone", "Brain", "GAT", "Heart", "Kidney", "Limb_Muscle", "Liver", "Lung", "Marrow", "MAT",
          "Pancreas", "SCAT", "Skin", "Small_Intestine", "Spleen", "WBC"),
  protein = c("kidney", "hippocampus", "heart", "fat", "cerebellum", "striatum", "spleen", "skeletalmuscle", "lung", "liver"),
  not_show = NULL
)

# change analytical design based on the data source (RNA/Protein)
choices_design <- list(
  rna = c("linear model using all age groups" = "linear",
          "two age groups comparison" = "two_group"),
  protein = c("two age groups comparison" = "two_group")
)

# Age group can be chosen in bulkRNA dataset
ages_rna <- c(1, 3, 6, 9, 12, 15, 18, 21, 24, 27)
