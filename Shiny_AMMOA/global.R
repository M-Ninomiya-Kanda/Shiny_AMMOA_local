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

library(clusterProfiler)
library(pathview)

library(tidyverse)

# load supporting file (mostly in-house functions)
source("shiny_multiomic_source.R")

# Load Necessary Data -----

# result of metabolomic estimation result from each tissue (CellMet 2025)
dat_metab <- readRDS("data/Metabolome/CellMet2025_metabolites_adjustedLFC_2026-09-13_exported.rds")
# compound IDs annotations
dat_compound_id <- read.csv("data/Metabolome/compound_mapping_res_2025-11-25.csv")

# load KEGG pathway IDs
kegg_ids <- readRDS("data/KEGG_pathway_ID_list_2026_09_15_exported.rds")
dat_pathway2compound <- readRDS("data/Metabolome/kegg_pathway2compound_list.RDS")
dat_pathway2name <- readRDS("data/Metabolome/kegg_pathway2name_list.RDS")

# Sample metadata
metadata_table_rna <- readRDS("data/Transcriptome/bulkTMS_sample_size_table_2026_09_15.rds")
metadata_table_prot <- readRDS("data/Proteome/Proteome_sample_size_table_2026_09_15.rds")
metadata_table_met <- readRDS("data/Metabolome/Metabolome_sample_size_table_2026_09_15.rds")

# Define Objects Required for Dynamic UI Change -----
# change tissue list in response to datasource (RNA/Protein)
choices_tissue <- list(
  rna = c("BAT", "Bone", "Brain", "GAT", "Heart", "Kidney", "Limb_Muscle", "Liver", "Lung", "Marrow", "MAT",
          "Pancreas", "SCAT", "Skin", "Small_Intestine", "Spleen", "WBC"),
  protein = c("kidney", "hippocampus", "heart", "fat", "cerebellum", "striatum", "spleen", "skeletalmuscle", "lung", "liver"),
  not_show = NULL
)

# change analytical design based on the EA algorithm
choices_design <- list(
  ora = c("across all ages" = "linear",
          "two age groups comparison" = "two_group"),
  gsea = c("across all ages" = "linear")
)

# Age group can be chosen in bulkRNA dataset
ages_rna <- c(1, 3, 6, 9, 12, 15, 18, 21, 24, 27)
