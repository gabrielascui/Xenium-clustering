#!/usr/bin/R

# ---
# Author: Francisco Emmanuel Castaneda-Castro
# Other contributors: 
    # Donaldo Sosa-Garcia
# Taking as a base Ciro's single-cell clustering pipeline (https://github.com/vijaybioinfo/clustering)
# Date: 2025-06-24
# Version 1.0.0
    # This code does the DGEA per resolution. 
# ---

## As a first step we need to open an interactive job hehe sorry.
#srun --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=300g --time=100:00:00 --pty bash -i
# module load R/4.3.3
# R

######################
# Clustering: Seurat: Spatial analysis #
######################
lib4.3path = "/home/fcastaneda/R/x86_64-pc-linux-gnu-library/4.3"
if(grepl("4.3", getRversion()) && file.exists(lib4.3path))
  .libPaths(new = lib4.3path)
options(future.globals.maxSize= 20480*1024^2)

library(optparse)
optlist <- list(
  optparse::make_option(
    opt_str = c("-y", "--yaml"), type = "character",
    help = "Configuration file: Instructions in YAML format."
  ),
  optparse::make_option(
    opt_str = c("-p", "--percent"), type = "numeric",
    help = "Percentage of variance."
  ),
  optparse::make_option(
    opt_str = c("-c", "--chosen_comp"), type = "numeric",
    help = "Chosen components."
  ),
    optparse::make_option(
    opt_str = c("-r", "--resolution"), type = "numeric", default = 50,
    help = "Resolution chosen."
  ),
  optparse::make_option(
    opt_str = c("-n", "--prefix"), type = "character",
    help = "Prefix for output. Name of the files."
  ),
  optparse::make_option(
    opt_str = c("-i", "--init_file"), type = "character",
    help = "Prefix for output. Name of the files."
  ),
  optparse::make_option(
    opt_str = c("--init"), type = "character",
    help = "Init object file from Snakemake."
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
library(BPCells)

sample_cells_per_cluster <- function(object, cluster_column, max_cells_per_cluster = 5000, seed = 1) {
  if(is.null(max_cells_per_cluster) || is.na(max_cells_per_cluster) || max_cells_per_cluster <= 0) {
    return(Cells(object))
  }
  set.seed(seed)
  cells_by_cluster <- split(rownames(object@meta.data), object@meta.data[, cluster_column])
  unname(unlist(lapply(cells_by_cluster, function(cluster_cells) {
    if(length(cluster_cells) > max_cells_per_cluster) {
      sample(cluster_cells, max_cells_per_cluster)
    } else {
      cluster_cells
    }
  })))
}

# opt parameters have the priority
if(interactive()){ # Example/manually
  opt$yaml = "config.yaml"
  opt$verbose <- TRUE
  opt$percent <- 50
  opt$resolution <- 0.8
  opt$chosen_comp <- 30
  opt$init_file <- ".object_init_mean0.01_pct50_pc30.rds"
  opt$prefix <- "init_mean0.01_pct50_pc30_res0.8"
}

config = yaml::read_yaml(opt$yaml)
cellranger_out= config$input_expression
pcts= opt$percent
pc= opt$chosen_comp
res <- opt$resolution
max_cells_per_cluster <- config$markers$max_cells_per_cluster
if(is.null(max_cells_per_cluster)) max_cells_per_cluster <- 5000
marker_seed <- config$markers$seed
if(is.null(marker_seed)) marker_seed <- 1
if(!is.null(opt$init) && is.null(opt$init_file)) opt$init_file <- opt$init

if(interactive()) outdir_sp= paste0(config$output_dir, "/", config$project_name)
outdir_sp = "."
cat("Working in:", getwd(), "\n")

if(opt$verbose) cat('Date and time:\n') ; st.time <- timestamp();
## Chekinf if parquet files are already in csv as required for Seurat
if(!dir.exists(outdir_sp)) dir.create(outdir_sp, recursive=TRUE)
setwd(outdir_sp)
# Load the init object --> Xenium

checking_init <- paste0("pct", pcts, "_pc", pc, ".rds")
if(!grepl(checking_init, opt$init_file)) stop("Init file doesn't mach with pct or pc given!")
xenium.obj <- readRDS(opt$init_file)

      for(resolution_cho in config$resolution) {
        comn_resolution<<-paste0(pcts, "pct", pc, "pc_res",resolution_cho)
        dir.create(comn_resolution, showWarnings = FALSE)
        cat(paste("Performing: ",comn_resolution, "\n"))
        cluster_value<- ifelse(grepl(config$norm, "sctransform"), paste0("SCT_snn_res.", resolution_cho), paste0("Xenium_snn_res.", resolution_cho))

        xenium.obj@meta.data[, cluster_value] <- factor(xenium.obj@meta.data[, cluster_value], levels=c(0:length(unique(xenium.obj@meta.data[, cluster_value]))))

        Idents(xenium.obj) <- cluster_value
        marker_obj <- xenium.obj
        if(grepl("MAST", config$markers$test, ignore.case = TRUE)) {
          marker_cells <- sample_cells_per_cluster(
            object = xenium.obj,
            cluster_column = cluster_value,
            max_cells_per_cluster = max_cells_per_cluster,
            seed = marker_seed
          )
          cat("MAST marker testing using ", length(marker_cells), " sampled cells; max ",
              max_cells_per_cluster, " per cluster.\n", sep = "")
          marker_obj <- subset(xenium.obj, cells = marker_cells)
          Idents(marker_obj) <- cluster_value
        }
           cat("Running DGEA with FindAllMarkers \n")

          markers_file_name <- paste0("/4.dgea_",config$markers$test, "_fc",config$markers$avg_logFC,"_padj", config$markers$p_val_adj)
          file_dgea <- paste0(comn_resolution, "/", markers_file_name, "_summary_stats_complete.csv")
          file_dgea_1 <- paste0(comn_resolution, "/.", markers_file_name, "_summary_stats.csv")
          
          if(!file.exists(file_dgea)){ 
            all_markers<-FindAllMarkers(object = marker_obj, 
              only.pos = FALSE,
              min.pct = 0.25,
              logfc.threshold = 0,
              min.diff.pct = 0.05,
              test.use = config$markers$test,
              return.thresh = 0.2,
              verbose = TRUE
            )

            cmarkers <- markers_summary(marktab = all_markers,
                  annot = marker_obj@meta.data,
                  datavis = GetAssayData(marker_obj),
                  cluster_column = cluster_value,
                  datatype = "SeuratNormalized",
                  verbose = opt$verbose)

            cmarkers_f1 <- cmarkers <- cmarkers[cmarkers$avg_log >= config$markers$avg_logFC & cmarkers$p_val_adj <= config$markers$p_val_adj, ]
              # tvar <- all_markers$avg_log >= 0.5 & all_markers$p_val_adj <= 0.05
            write.csv(cmarkers_f1, file_dgea_1)

            void <- marker_report(
              markers_df = cmarkers_f1,
              file = "_summary_stats.csv",
              verbose = TRUE, return_plot = TRUE
            )

            pdf(paste0(comn_resolution, markers_file_name, "deltaSlope.pdf"), height=14, width=20)
            print(void[[1]])
            dev.off()

            pdf(paste0(comn_resolution, markers_file_name, "delta_percentage.pdf"), height=14, width=20)
            print(void[[2]])
            dev.off()
            
            genes_to_sum<-rownames(marker_obj)[!rownames(marker_obj) %in% cmarkers$gene]

            to_join_stats_gene<-stats_summary_table(
                mat = GetAssayData(marker_obj),
                groups = make_list(x = marker_obj@meta.data, colname = cluster_value, grouping = TRUE),
                rnames = genes_to_sum,
                datatype = "SeuratNormalized",
                verbose = TRUE
            )
            
            cmarkers_complete <- bind_rows(cmarkers, to_join_stats_gene)

            write.csv(cmarkers_complete, file_dgea)
          } else {cat ("Taking DGEA from file \n"); cmarkers_f1 <- read.csv(file_dgea_1)} 

          aver<-cmarkers_f1 %>% group_by(cluster) %>% top_n(8, Dpct)
          
          genes_plot_out<-paste0(comn_resolution, "/5.genes_examples")
          dir.create(genes_plot_out)

          aver <- aver %>% select(cluster, gene)
          custom_genes <- data.frame(
            cluster = "_custom", 
            gene = config$markers$to_plot
          )

          aver <- rbind(aver, custom_genes)
          # aver %>% print(n=100)

          lapply(unique(aver$cluster), function(x){
            # x<-1
            feat_clust<- aver %>% filter(cluster == x) %>% pull(gene)
            # feat_clust <- c(feat_clust)
            pdf(paste0(genes_plot_out, "/cluster", x, ".pdf"))
            repet_plot<-c(0:ceiling(length(feat_clust)/4))  ### Check Donas genes modi
            lapply(repet_plot, function(nplot){ 
              nplot<-nplot*4
              try(print(FeaturePlot(xenium.obj,features = feat_clust[nplot:(nplot+3)], max.cutoff="q10", order = TRUE)))
              })
            dev.off()
            }) 
          rm(marker_obj)
      }

        rm(xenium.obj); rm(cmarkers_f1); rm(aver)

      
if(opt$verbose) cat("Saving output file for snakemake \n")
done_file <- paste0(".", gsub("init", "markers", opt$prefix), ".txt")
writeLines("markers done", done_file)

if(opt$verbose){
  cat('\n\n*******************************************************************\n')
  cat('Starting time:\n'); cat(st.time, '\n')
  cat('Finishing time:\n'); timestamp()
  cat('*******************************************************************\n')
  cat('SESSION INFO:\n'); print(sessionInfo()); cat("\n")
  cat('Pipeline finished successfully\n')
}
  
