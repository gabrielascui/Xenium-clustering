#!/usr/bin/R

# ---
# Author: Francisco Emmanuel Castaneda-Castro
# Contributors: 
    # Donaldo Sosa-Garcia
# Taking as a base Ciro's single-cell clustering pipeline (https://github.com/vijaybioinfo/clustering)
# Date: 2025-06-24
# Version 1.0.0
    # Seurat object creation, first QC's, metadata filtering, normalization and high variable genes (HVG slection), PC and UMAP
    #Version 2.0.0:
        # Performing the analysis for multiple slides based on Seurat reocmendations. Using the merge fucntion. 
        # Plots improvements and plot for each slide independent
# ---

######################
# Clustering: Seurat: Spatial analysis #
######################
# plan("multiprocess") ## Use multicore
lib4.3path = "/home/fcastaneda/R/x86_64-pc-linux-gnu-library/4.3"
if(grepl("4.3", getRversion()) && file.exists(lib4.3path))
  .libPaths(new = lib4.3path)
library(future)
options(future.globals.maxSize= 26312*1024^2)
# library(future)
# plan("multisession", workers=availableCores())

optlist <- list(
  optparse::make_option(
    opt_str = c("-y", "--yaml"), type = "character",
    help = "Configuration file: Instructions in YAML format."
  ),
    optparse::make_option(
    opt_str = c("-v", "--verbose"), default = TRUE,
    help = "Verbose: Show progress."
  )
)

# Getting arguments from command line
opt <- optparse::parse_args(optparse::OptionParser(option_list = optlist))

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if(length(script_file) && !is.na(script_file)) dirname(normalizePath(script_file)) else "R"
source(file.path(script_dir, "functions.R"))

library(ggplot2)
library(reshape2)
library(dplyr)
library(Seurat)
library(parallel)
library(Seurat)
library(plotly)
library(stringr)
library(ggrepel)
library(optparse)
library(cowplot)
library(BPCells)

cat("New version \n")
# opt parameters have the priority
if(interactive()){ # Example/manuallyxw
  opt$yaml = "config.yaml"
  opt$verbose <- TRUE
}

config = yaml::read_yaml(opt$yaml)
input_expression= config$input_expression
if(is.character(input_expression) && length(input_expression) == 1) input_expression <- as.list(input_expression)
cellranger_out = config$input_expression
if(is.character(cellranger_out) && length(cellranger_out) == 1) cellranger_out <- as.list(cellranger_out)
pcts= opt$percent
pc= opt$chosen_comp

outdir_sp = "."
if(interactive()) outdir_sp= paste0(config$output_dir, "/", config$project_name)

cat("Working in:", getwd(), "\n")


if(opt$verbose) cat('Date and time:\n') ; st.time <- timestamp();
## Chekinf if parquet files are already in csv as required for Seurat
if(!dir.exists(outdir_sp)) dir.create(outdir_sp, recursive=TRUE)
setwd(outdir_sp)
# Load the Xenium data

#################################################
################# Preprocessing #################
#################################################
outdir_final <- paste0(config$output_dir, "/", config$project_name)

attach_bpcells_counts <- function(xenium_obj, bp_dir) {
  if (!dir.exists(bp_dir)) {
    cat("Writing BPCells matrix:", bp_dir, "\n")
    write_matrix_dir(mat = xenium_obj[["Xenium"]]$counts, dir = bp_dir)
  } else {
    cat("Using existing BPCells matrix:", bp_dir, "\n")
  }
  counts.mat <- open_matrix_dir(dir = bp_dir)
  assay_to_upd <- CreateAssay5Object(counts = counts.mat, key = xenium_obj@assays$Xenium@key)
  xenium_obj@assays$Xenium@layers <- list()
  xenium_obj@assays$Xenium <- assay_to_upd
  xenium_obj[["Xenium"]]$counts <- counts.mat
  xenium_obj
}

load_slide_as_bpcells <- function(cellranger_path, slide_index) {
  slide_name <- paste0("S", slide_index)
  cat("Starting with ", cellranger_path, "\n")
  xen_l <- LoadXenium(cellranger_path)
  xen_l <- RenameCells(xen_l, new.names = paste0(slide_name, "_", Cells(xen_l)))
  xen_l$origlib <- slide_name
  names(xen_l@images) <- paste0("fov.", slide_index)
  bp_dir <- file.path(outdir_final, paste0("BPcells_matrix_slide", slide_index))
  xen_l <- attach_bpcells_counts(xen_l, bp_dir)
  saveRDS(xen_l, file.path(outdir_final, paste0(".object_input_slide", slide_index, ".rds")))
  cat(slide_name, " done\n")
  xen_l
}

input_object_file <- file.path(outdir_final, ".object_input_data.rds")
if(file.exists(input_object_file)){
    cat("Reading saved matrix file in .object_input_data.rds \n")
    xenium.obj <- readRDS(input_object_file)
    if (file.normalize(input_object_file, mustWork = FALSE) != file.normalize(".object_input_data.rds", mustWork = FALSE)) {
      file.copy(input_object_file, ".object_input_data.rds", overwrite = TRUE)
    }
    cat("Finished \n")
} else {
  if(length(input_expression) ==  1 & all(grepl(".rds", unlist(input_expression)))){ 
    cat("Reading input expression matrix: If you provide to the pipeline an rds file with multiple slides I'm assuming they are already merged. \n ")
    if(file.exists(unlist(input_expression))){  
      xenium.obj <- readRDS(unlist(input_expression))
      xenium.obj <- attach_bpcells_counts(xenium.obj, file.path(outdir_final, "BPcells_matrix_slide1"))
    } else { 
    stop("Error file doesn't exists") 
    }
  } else {
    if (length(cellranger_out) == 1) {
      cat("------- Reading object -------- \n")
      xenium.obj <- load_slide_as_bpcells(cellranger_out[[1]], 1)
    } else {
      cat('  -------  Reading xenium datasets: more than one Xenium dataset found  -------  \n')
      cat(paste0('Xenium datasets found: ', length(cellranger_out), '\n'))

      xenium.obj_list <- lapply(seq_along(cellranger_out), function(x) {
        load_slide_as_bpcells(cellranger_out[[x]], x)
      })

      cat('  -------  Merging seurat objects  -------  \n')
      xenium.obj <- Reduce(function(x, y) merge(x, y), xenium.obj_list)
      rm(xenium.obj_list)
    }  
  }

  perc_genes_keep <- config$filtering$perc_genes_keep
  if(is.null(perc_genes_keep)) perc_genes_keep <- 100
  if(perc_genes_keep != 100){ 
    # config$filtering$perc_genes_keep <- 20

      if(opt$verbose) cat('\n ######### Subseting features accordingly to config file ...\n Keeping only ', config$filtering$perc_genes_keep, '% top of the genes. Selecting the genes based on their variance'); timestamp()
      
      if(opt$verbose) cat('\n@@@@@@@@@ Selecting features (FindVariableFeatures) ...\n'); timestamp()
      hvg_final <- rownames(xenium.obj)
      xenium.obj <- FindVariableFeatures(
        object = xenium.obj,
        selection.method = config$variable_features$method,
        nfeatures = config$variable_features$nfeatures,
        mean.cutoff = config$variable_features$mean.cutoff,
        dispersion.cutoff =  config$variable_features$dispersion.cutoff,
        verbose = opt$verbose
      )
      hvg_df <- HVFInfo(xenium.obj, method = "vst")
      disp_n <- grep("standardized", colnames(hvg_df), value=TRUE); hvg_df <- hvg_df[order(-hvg_df[, disp_n]), ] # order
      # mean_pct_filters <- TRUE

      if(opt$verbose) cat("Taking", perc_genes_keep, "%\n")
      hvg_df_keep<-hvg_df
      hvg_df_keep$cumulative = cumsum(hvg_df_keep[, disp_n])
      hvg_df_keep$cumulative_pct = round(hvg_df_keep$cumulative / sum(hvg_df_keep[, disp_n], na.rm = TRUE) * 100, 2)
      passed <- hvg_df_keep$cumulative_pct <= perc_genes_keep
      if(opt$verbose) cat("Number of features:", sum(passed), "of", nrow(hvg_df_keep), "\n")
      hvg_df_keep$passed<- hvg_df_keep$cumulative_pct <= perc_genes_keep & (rownames(hvg_df_keep) %in% hvg_final) ## The last condition could be  when variable "passed" is definded
      hvg_df_keep <- hvg_df_keep[passed, ]; #log_history <- paste0(log_history, "Pct")
      hvg_df_keep <- rownames(hvg_df_keep)

      not_selected_genes <- rownames(xenium.obj)[!rownames(xenium.obj) %in% hvg_df_keep]
      
      recover_genes <- rowSums(xenium.obj@assays$Xenium$counts[not_selected_genes,colSums(xenium.obj@assays$Xenium$counts[hvg_df_keep,]) == 0])

      biggest_genes_recovered <- sort(recover_genes, decreasing=TRUE)[1:200]
      hvg_df_keep2 <- unique(c(hvg_df_keep, names(biggest_genes_recovered)))
      if(opt$verbose) cat("Keeping", length(hvg_df_keep2), "of", nrow(hvg_df), " features:\n")

      keep_cells <- colnames(xenium.obj)[colSums(xenium.obj@assays$Xenium$counts[hvg_df_keep2,]) > 0]
      if(opt$verbose) cat("Keeping", length(keep_cells), "of", ncol(xenium.obj), " cells:\n")

      # biggest_genes <- sort(rowSums(xenium.obj@assays$Xenium$counts[hvg_df_keep,]), decreasing=FALSE)[1:100]
      # hvg_df_keep2 <- unique(c(hvg_df_keep, names(biggest_genes)))
      xenium.obj <- subset(xenium.obj, cells=keep_cells, features = hvg_df_keep2)
  }

      cat("Saving matrix file in .object_input_data.rds \n")
      saveRDS(xenium.obj, ".object_input_data.rds")
      cat("Finished \n")

} 

if(opt$verbose){
  cat('\n\n*******************************************************************\n')
  cat('Starting time:\n'); cat(st.time, '\n')
  cat('Finishing time:\n'); timestamp()
  cat('*******************************************************************\n')
  cat('SESSION INFO:\n'); print(sessionInfo()); cat("\n")
  cat('Pipeline finished successfully\n')
}
  
