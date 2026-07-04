#!/usr/bin/R

######################
# Clustering: Seurat: Spatial analysis
######################

options(future.globals.maxSize= 26312*1024^2)
lib4.3path = "/home/fcastaneda/R/x86_64-pc-linux-gnu-library/4.3"
if(grepl("4.3", getRversion()) && file.exists(lib4.3path))
  .libPaths(new = lib4.3path)

ElbowPlot2 <- function (object, ndims = 20, reduction = "pca"){
  data.use <- Stdev(object = object, reduction = reduction)
  if (length(x = data.use) == 0) {
      stop(paste("No standard deviation info stored for", reduction))
  }
  if (ndims > length(x = data.use)) {
      warning("The object only has information for ", length(x = data.use), " reductions")
      ndims <- length(x = data.use)
  }
  stdev <- "Standard Deviation"
  plot <- ggplot(data = data.frame(dims = 1:ndims, stdev = data.use[1:ndims])) +
      geom_point(mapping = aes_string(x = "dims", y = "stdev"), size = 2) +
      labs(x = gsub(pattern = "_$", replacement = "", x = Key(object = object[[reduction]])),
          y = stdev) + theme_cowplot()
  return(plot)
}

optlist <- list(
  optparse::make_option(c("-y", "--yaml"), type = "character", help = "Configuration file: Instructions in YAML format."),
  optparse::make_option(c("-s", "--stage"), type = "character", default = "cluster",
    help = "Pipeline stage: normalize_base, hvg, pca, or cluster."),
  optparse::make_option(c("-p", "--percent"), type = "numeric", help = "Percentage of variance."),
  optparse::make_option(c("-n", "--n_comp"), type = "numeric", default = 50, help = "Total number of components to explore."),
  optparse::make_option(c("-c", "--chosen_comp"), type = "numeric", help = "Chosen components."),
  optparse::make_option(c("-r", "--prefix"), type = "character", help = "Prefix for output. Name of the files."),
  optparse::make_option(c("-v", "--verbose"), default = TRUE, help = "Verbose: Show progress.")
)

opt <- optparse::parse_args(optparse::OptionParser(option_list = optlist))

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if(length(script_file) && !is.na(script_file)) dirname(normalizePath(script_file)) else "R"
source(file.path(script_dir, "functions.R"))

library(ggplot2)
library(reshape2)
library(dplyr)
library(Seurat)
library(reticulate)
Sys.setenv(RETICULATE_PYTHON = "~/miniconda3/envs/clustering_new/bin/python")
library(reticulate)
py_module_available(module = "leidenalg")
library(plotly)
library(stringr)
library(ggrepel)
library(optparse)
library(cowplot)
library(BPCells)

# RunPCA can coerce BPCells matrices through Matrix limits. This keeps variance
# calculation on the IterableMatrix path when possible.
fixed_PrepDR5 <- function(object, features = NULL, layer = "scale.data", verbose = TRUE) {
  layer <- layer[1L]
  olayer <- layer
  layer <- SeuratObject::Layers(object = object, search = layer)
  if (is.null(layer)) {
    abort(paste0("No layer matching pattern '", olayer, "' not found. Please run ScaleData and retry"))
  }
  data.use <- SeuratObject::LayerData(object = object, layer = layer)
  features <- features %||% VariableFeatures(object = object)
  if (!length(x = features)) {
    stop("No variable features, run FindVariableFeatures() or provide a vector of features", call. = FALSE)
  }
  if (is(data.use, "IterableMatrix")) {
    features.var <- BPCells::matrix_stats(matrix=data.use, row_stats="variance")$row_stats["variance",]
  } else {
    features.var <- apply(X = data.use, MARGIN = 1L, FUN = var)
  }
  features.keep <- features[features.var > 0]
  if (!length(x = features.keep)) {
    stop("None of the requested features have any variance", call. = FALSE)
  } else if (length(x = features.keep) < length(x = features)) {
    exclude <- setdiff(x = features, y = features.keep)
    if (isTRUE(x = verbose)) {
      warning("The following ", length(x = exclude),
        " features requested have zero variance; running reduction without them: ",
        paste(exclude, collapse = ", "), call. = FALSE, immediate. = TRUE)
    }
  }
  features <- features.keep
  features <- features[!is.na(x = features)]
  features.use <- features[features %in% rownames(data.use)]
  if(!isTRUE(all.equal(features, features.use))) {
    missing_features <- setdiff(features, features.use)
    if(length(missing_features) > 0) {
      warning(paste("The following features were not available: ",
        paste(missing_features, collapse = ", "), ".", sep = ""), immediate. = TRUE)
    }
  }
  data.use[features.use, ]
}

config = yaml::read_yaml(opt$yaml)
outdir_sp = "."
if(interactive()) outdir_sp = paste0(config$output_dir, "/", config$project_name)
if(!dir.exists(outdir_sp)) dir.create(outdir_sp, recursive=TRUE)
setwd(outdir_sp)

cat("Working in:", getwd(), "\n")
if(opt$verbose) cat("Date and time:\n")
st.time <- timestamp()

object_file <- function(prefix) paste0(".object_", prefix, ".rds")
base_file <- ".object_normalized_base.rds"
input_file <- ".object_input_data.rds"

norm_prefix_from_pca <- function(prefix) sub("^pca", "norm", prefix)
pca_prefix_from_init <- function(prefix) sub("_pc[^_]+$", "", sub("^init", "pca", prefix))

read_object <- function(path) {
  if(!file.exists(path)) stop("Required input object not found: ", path)
  readRDS(path)
}

run_qc_and_metadata <- function(xenium.obj) {
  cat("------- Overall QCs -------- \n")
  dir.create("Section_qc", showWarnings = FALSE)

  pdf("Section_qc/TissueSection_Prefilter.pdf", width=10, height=5)
  lapply(names(xenium.obj@images), function(slide_fov){
    ImageFeaturePlot(xenium.obj, features = c("nCount_Xenium"), fov = slide_fov,
      border.size=NA, max.cutoff=400, axes=TRUE)
  })
  dev.off()

  meta_data <- xenium.obj@meta.data
  to_text <- meta_data %>% group_by(origlib) %>% summarize(size=n())

  pdf(paste0("Section_qc/QCs_prefiltered_", ncol(xenium.obj), "_cells.pdf"), width=5, height=7)
  vlns <- lapply(c("nCount_Xenium", "nFeature_Xenium"), function(x){
    qc_violin_spatial(dat = meta_data, yax = x, xax = "origlib", to_text=to_text)
  })
  print(cowplot::plot_grid(plotlist = vlns, ncol = 1))
  dev.off()

  pdf(paste0("Section_qc/QCs_prefiltered_ControlProbes_", ncol(xenium.obj), "_cells.pdf"), width=5, height=19)
  vlns <- lapply(c("nCount_ControlProbe", "nFeature_ControlProbe","nCount_ControlCodeword",
    "nFeature_ControlCodeword", "nCount_BlankCodeword", "nFeature_BlankCodeword"), function(x){
    qc_violin_spatial(dat = meta_data, yax = x, xax = "origlib", to_text=to_text)
  })
  print(cowplot::plot_grid(plotlist = vlns, ncol = 1))
  dev.off()
  rm(vlns)

  pdf("Section_qc/UMI_per_gene_distribution.pdf", width=10, height=5)
  dist <- lapply(names(xenium.obj@assays$Xenium@layers), function(slide_fov) {
    umi_per_genes <- rowSums(xenium.obj@assays$Xenium@layers[[slide_fov]]) %>% melt() %>% arrange(value)
    write.csv(umi_per_genes, paste0("Section_qc/UMIs_per_gene_in_slide", gsub("counts\\.", "", slide_fov), ".csv"))
    slide_name <- gsub("counts.", "", slide_fov)
    ggplot(umi_per_genes, aes(x = log10(value))) + geom_density(color=4) +
      labs(x= "log10(UMI count per gene)") + ggtitle(paste0("slide", slide_name)) + theme_classic()
  })
  print(dist); dev.off(); rm(dist)

  pdf("Section_qc/Distribution-UMI_per_cell.pdf", width=10, height=5)
  dist <- lapply(unique(xenium.obj$origlib), function(lib_slide){
    md1 <- xenium.obj@meta.data[xenium.obj$origlib == lib_slide, ]
    umi_cells <- melt(table(md1$nCount_Xenium))
    xliC <- round(mean(md1$nCount_Xenium))
    ggplot(umi_cells, aes(y=value, x=Var1)) +
      geom_bar(color="blue", fill="blue", stat="identity") + labs(y= "Number of cells", x="#UMI") +
      theme_classic() + geom_vline(xintercept = xliC, linetype="dotted", color = "red", size=1.5) +
      labs(subtitle = paste("mean nCount: ", xliC)) + ggtitle(paste0("Slide ", lib_slide))
  })
  print(dist); dev.off(); rm(dist)

  pdf("Section_qc/Distribution-genes_per_cell.pdf", width=10, height=5)
  dist <- lapply(unique(xenium.obj$origlib), function(lib_slide){
    md1 <- xenium.obj@meta.data[xenium.obj$origlib == lib_slide, ]
    genes_cells <- melt(table(md1$nFeature_Xenium))
    xliF <- round(mean(md1$nFeature_Xenium))
    xliC <- sum(md1$nFeature_Xenium > xliF)
    ggplot(genes_cells, aes(y=value, x=Var1)) +
      geom_bar(color="blue", fill="blue", stat="identity") +
      labs(y= "Number of cells", x="# genes with > 3 UMI") + theme_classic() +
      geom_vline(xintercept = xliF, linetype="dotted",  color = "red", size=1.5) +
      labs(subtitle = paste("mean nFeatures: ", xliF, ". \nCells with more than mean: ", xliC)) +
      ggtitle(paste0("Slide ", lib_slide))
  })
  print(dist); dev.off(); rm(dist)

  meta_data <- xenium.obj@meta.data

  if(!is.null(config$metadata)) {
    if(opt$verbose) cat(" ------- Adding metadata -------  \n")
    addmetadataf <- if(length(config$metadata) > 0) config$metadata else "no_file"
    addmetadataf <- path.expand(addmetadataf)
    tvar <- file.exists(addmetadataf)
    if(any(tvar) & length(tvar) > 1){
      addmetadataf <- addmetadataf[tvar]
      slide_cell_md <- lapply(seq_along(addmetadataf), function(x) {
        cell_md <- addmetadataf[[x]]
        addannot <- if(grepl(".csv", cell_md)) remove.factors(readfile(cell_md, row.names=1)) else remove.factors(readfile(cell_md))
        cells_per_slide <- table(meta_data[gsub("S.*_", "", rownames(meta_data)) %in% rownames(addannot),"origlib"])/table(meta_data$origlib)*100
        if(!any(cells_per_slide > 30)) warning("A really low proportion of cells match with a slide.")
        slide_name <- names(cells_per_slide)[rev(order(cells_per_slide))[1]]
        rownames(addannot) <- paste0(slide_name, "_", rownames(addannot))
        addannot$cell_id <- paste0(slide_name, "_", addannot$cell_id)
        addannot$orig.TMA_core_slide <- paste0(slide_name, "_", addannot$orig.TMA_core)
        addannot
      })
      slide_cell_md2 <- Reduce(rbind, slide_cell_md)
      meta_data <- meta_data[rownames(meta_data) %in% rownames(slide_cell_md2),]
      meta_data <- joindf(x = meta_data, y = slide_cell_md2)
    } else if (any(tvar) & length(tvar) == 1) {
      addmetadataf <- addmetadataf[tvar]
      addannot <- if(grepl(".csv", addmetadataf)) remove.factors(readfile(addmetadataf, row.names=1)) else remove.factors(readfile(addmetadataf))
      meta_data <- meta_data[rownames(meta_data) %in% rownames(addannot),]
      meta_data <- joindf(x = meta_data, y = addannot)
    }
  } else cat("No extra metadata given \n")

  if(!is.null(config$metadata_donor)){
    if(opt$verbose) cat(" ------- Adding donor metadata -------  \n")
    tvar <- as.character(meta_data$orig.TMA_core)
    maxln <- max(sapply(tvar, length))
    meta_donor <- data.table::rbindlist(lapply(tvar, function(x){
      if(length(x) < maxln && length(x) == 1) x <- rep(x, length.out = maxln)
      if(length(x) < maxln) x <- rep(paste0(x, collapse = "-"), length.out = maxln)
      as.data.frame(t(x))
    }))
    meta_donor <- remove.factors(data.frame(meta_donor, row.names = rownames(meta_data)))
    rnname <- ifelse(grepl("~", config$metadata_donor), gsub(".*~", "", config$metadata_donor), 1)
    config$metadata_donor <- gsub("~.*", "", config$metadata_donor)

    if(length(config$metadata_donor) > 1) {
      mdonor_extra2 <- lapply(seq_along(config$metadata_donor), function(x) {
        md_donor_s <- config$metadata_donor[[x]]
        tmp <- read.csv(md_donor_s, stringsAsFactors = FALSE, row.names = rnname[[x]])
        rownames(tmp) <- paste0("S", x, "_", rownames(tmp))
        tmp[tmp == ""] <- NA
        md_donor_merged <- tmp[meta_donor$V1[meta_donor$V1 %in% rownames(tmp)], , drop = FALSE]
        rownames(md_donor_merged) <- rownames(meta_donor[meta_donor$V1 %in% grep(paste0("^S", x), meta_donor$V1, value = TRUE), , drop =FALSE])
        md_donor_merged[] <- lapply(md_donor_merged, function(col) paste0("S", x, "_", col))
        md_donor_merged$cells <- rownames(md_donor_merged)
        md_donor_merged
      })
      mdonor_extra <- Reduce(function(x, y) {
        to_joincols <- intersect(names(x), names(y))
        full_join(x, y, by = to_joincols)}, mdonor_extra2)
      rownames(mdonor_extra) <- mdonor_extra$cells
      meta_donor <- joindf(data.frame(meta_donor), mdonor_extra)
    } else {
      tmp <- read.csv(config$metadata_donor, stringsAsFactors = FALSE, row.names = rnname)
      tmp[tmp == ""] <- NA
      mdonor_extra <- tmp[meta_donor$V1, , drop = FALSE]
      rownames(mdonor_extra) <- rownames(meta_donor)
      meta_donor <- joindf(data.frame(meta_donor), mdonor_extra)
    }
    colnames(meta_donor) <- paste0("orig.", colnames(meta_donor))
    meta_data <- joindf(meta_data, meta_donor)
  } else cat("No extra metadata donor given \n")

  if(!is.null(config$filtering$subset)) {
    cat("  -------  Filters present  -------  \n")
    filtereddata <- filters_complex(
      mdata = meta_data,
      filters = lapply(names(config$filtering$subset), function(x) c(x, config$filtering$subset[[x]]) ),
      verbose = opt$verbose
    )
    cat("Preserving:", nrow(filtereddata[[1]]), "/", nrow(meta_data), "samples/cells\n")
    meta_data <- filtereddata[[1]]
    xenium.obj <- subset(xenium.obj, cells = rownames(meta_data))
    meta_data <- meta_data %>% select(!matches("seurat_clus|^Xenium_snn_res.|^SCT_snn_res.|umap_1|umap_2"))
    xenium.obj <- AddMetaData(object = xenium.obj, metadata = meta_data)
  } else cat("No filters present ")

  meta_data <- xenium.obj@meta.data
  to_text <- meta_data %>% group_by(origlib) %>% summarize(size=n())

  pdf(paste0("Section_qc/QCs_postfiltered_", ncol(xenium.obj), "_cells.pdf"), width=5, height=7)
  vlns <- lapply(c("nCount_Xenium", "nFeature_Xenium"), function(x){
    qc_violin_spatial(dat = meta_data, yax = x, xax = "origlib", to_text=to_text)
  })
  print(cowplot::plot_grid(plotlist = vlns, ncol = 1))
  dev.off()

  pdf(paste0("Section_qc/QCs_postfiltered_ControlProbes_", ncol(xenium.obj), "_cells.pdf"), width=5, height=19)
  vlns <- lapply(c("nCount_ControlProbe", "nFeature_ControlProbe","nCount_ControlCodeword",
    "nFeature_ControlCodeword", "nCount_BlankCodeword", "nFeature_BlankCodeword"), function(x){
    qc_violin_spatial(dat = meta_data, yax = x, xax = "origlib", to_text=to_text)
  })
  print(cowplot::plot_grid(plotlist = vlns, ncol = 1))
  dev.off()

  pdf("Section_qc/TissueSection_QCs_Postfiltered.pdf", width=10, height=5)
  for (fov in names(xenium.obj@images)){
    print(ImageFeaturePlot(xenium.obj, features = c("nCount_Xenium"), fov = fov, border.size=NA, max.cutoff=400))
    print(ImageFeaturePlot(xenium.obj, fov = fov, features = c("nFeature_Xenium"), border.size=NA))
  }
  dev.off()

  xenium.obj
}

run_normalize_base <- function() {
  xenium.obj <- read_object(input_file)
  xenium.obj <- run_qc_and_metadata(xenium.obj)

  hvg_final <- rownames(xenium.obj)
  tvar <- config$variable_features$file
  if(!is.null(tvar) && tvar != "no_file"){
    if(file.exists(tvar)){
      file_con <- file(description = tvar, open = "r")
      these_feats <- readLines(con = file_con)
      close(file_con)
      hvg_final <- hvg_final[hvg_final %in% these_feats]
    } else {
      stop("Error: file config$variable_features$file does not exist")
    }
  }

  if(!grepl(config$norm, "sctransform")){
    if(opt$verbose) cat("\n@@@@@@@@@ Normalizing (NormalizeData) ...\n")
    xenium.obj <- NormalizeData(
      object = xenium.obj,
      normalization.method = config$norm,
      scale.factor = 10000,
      verbose = opt$verbose
    )
    if(opt$verbose) cat("\n@@@@@@@@@ Selecting features (FindVariableFeatures) ...\n")
    xenium.obj <- FindVariableFeatures(
      object = xenium.obj,
      selection.method = config$variable_features$method,
      nfeatures = config$variable_features$nfeatures,
      mean.cutoff = config$variable_features$mean.cutoff,
      dispersion.cutoff =  config$variable_features$dispersion.cutoff,
      verbose = opt$verbose
    )
    xenium.obj@misc$hvg_final <- hvg_final
  } else {
    if(opt$verbose) cat("\n@@@@@@@@@ SCTransform - no regression\n")
    xenium.obj <- SCTransform(
      object = xenium.obj,
      variable.features.n = config$variable_features$nfeatures,
      vars.to.regress = NULL,
      conserve.memory = TRUE,
      return.only.var.genes = TRUE,
      verbose = opt$verbose
    )
  }
  saveRDS(object = xenium.obj, file = base_file)
}

run_hvg <- function() {
  xenium.obj <- read_object(base_file)

  if(!grepl(config$norm, "sctransform")){
    hvg_df <- HVFInfo(xenium.obj, method = "vst")
    disp_n <- grep("standardized", colnames(hvg_df), value=TRUE)
    hvg_df <- hvg_df[order(-hvg_df[, disp_n]), ]
    hvg_final <- xenium.obj@misc$hvg_final
    if(is.null(hvg_final)) hvg_final <- rownames(xenium.obj)
    mean_pct_filters <- !isTRUE(config$variable_features$file_only)

    if(mean_pct_filters){
      clustering_pct <- opt$percent
      if(opt$verbose) cat("Taking", clustering_pct, "%\n")
      hvg_df_keep <- hvg_df
      hvg_df_keep$cumulative = cumsum(hvg_df_keep[, disp_n])
      hvg_df_keep$cumulative_pct = round(hvg_df_keep$cumulative / sum(hvg_df_keep[, disp_n], na.rm = TRUE) * 100, 2)
      passed <- hvg_df_keep$cumulative_pct <= clustering_pct
      hvg_df_keep$passed <- passed & (rownames(hvg_df_keep) %in% hvg_final)

      pdf(paste0(outdir_sp, "/pct", clustering_pct, "_cumulativeVariance.pdf"))
      p <- ggplot(hvg_df_keep, aes(y=cumulative_pct, x=1:length(cumulative_pct), color=mean)) + geom_point() +
        scale_colour_gradient2(low="#ffdf32", mid="#ff9900", high="#ff5a00",
          midpoint=quantile(hvg_df_keep$mean, 0.95)/2, limits=c(0,quantile(hvg_df_keep$mean, 0.95))) +
        theme_classic() + geom_hline(yintercept = clustering_pct, linetype="dotted", col = "red", size = 1.3) +
        labs(title="Cumulative variance of Highly Variable Genes",
          x="Number of genes, order: variance.standardized",
          y="Cumulative percentage: variance standardized")
      print(p)
      dev.off()

      hvg_df_keep <- hvg_df_keep[passed, ]
      hvg_final2 <- intersect(hvg_final, rownames(hvg_df_keep))
      VariableFeatures(object = xenium.obj, assay = NULL) <- hvg_final2

      pdf(paste0(outdir_sp, "/pct", clustering_pct, "_VariableFeatures.pdf"))
      plot1 <- VariableFeaturePlotp_bicolor(object = xenium.obj, assay = "Xenium", variable_features = NULL)
      print(plot1)
      dev.off()
    }

    if(opt$verbose) cat("\n@@@@@@@@@ Scaling data...\n")
    xenium.obj <- ScaleData(
      object = xenium.obj,
      vars.to.regress = config$regress_var,
      block.size = 2000,
      verbose = opt$verbose
    )
  } else {
    if(opt$verbose) cat("Using SCTransform variable features; percent loop is ignored by SCT.\n")
  }

  saveRDS(object = xenium.obj, file = object_file(opt$prefix))
}

run_pca <- function() {
  input_norm <- object_file(norm_prefix_from_pca(opt$prefix))
  xenium.obj <- read_object(input_norm)

  assignInNamespace("PrepDR5", fixed_PrepDR5, "Seurat")
  if(opt$verbose){
    cat("\n@@@@@@@@@ Linear dimensional reduction\n")
    cat("Computing:", casefold(config$dim_reduction$base$type, upper = TRUE), "\n")
    cat("Components:", opt$n_comp, "\n")
  }

  xenium.obj <- RunPCA(xenium.obj,
    npcs = opt$n_comp,
    features = VariableFeatures(xenium.obj),
    nfeatures.print = 15,
    verbose = opt$verbose)

  if (grepl("harmony", config$dim_reduction$base$type, ignore.case = TRUE)) {
    xenium.obj <- harmony::RunHarmony(
      object = xenium.obj,
      group.by.vars = config$dim_reduction$base$batch,
      dims.use = 1:opt$n_comp,
      verbose = opt$verbose
    )
  }

  saveRDS(object = xenium.obj, file = object_file(opt$prefix))
}

run_cluster <- function() {
  input_pca <- object_file(pca_prefix_from_init(opt$prefix))
  xenium.obj <- read_object(input_pca)
  clustering_pct <- opt$percent
  clustering_per_pcs <- opt$chosen_comp
  clustering_per_pcs_num <- ifelse(clustering_per_pcs < 10, paste0("0", clustering_per_pcs), clustering_per_pcs)

  pdf(paste0("pct", clustering_pct, "_Elbow_", clustering_per_pcs, "pcs.pdf"), width=10, height=5)
  print(ElbowPlot2(xenium.obj, ndims = opt$n_comp, reduction = "pca") +
    geom_vline(xintercept = clustering_per_pcs, linetype="dotted", color = "red", size=1.5))
  dev.off()

  n_neis <- if(!is.null(config$dim_reduction$umap)) unique(config$dim_reduction$umap$n.neighbors)

  xenium.obj <- RunUMAP(xenium.obj,
    reduction = config$dim_reduction$base$type,
    dims = 1:clustering_per_pcs,
    n.neighbors = n_neis,
    min.dist = config$dim_reduction$umap$min.dist,
    verbose = opt$verbose)

  if(opt$verbose) cat("Resolutions:", paste0(config$resolution, collapse = ", "), "\n")

  xenium.obj <- FindNeighbors(xenium.obj,
    reduction = config$dim_reduction$base$type,
    dims = 1:clustering_per_pcs,
    verbose = opt$verbose)

  xenium.obj <- FindClusters(xenium.obj,
    resolution = config$resolution,
    verbose = opt$verbose)

  global_resolution <- paste0(ifelse(clustering_pct < 10, paste0("0", clustering_pct), clustering_pct),
    "pct", clustering_per_pcs_num, "pc_qc")
  dir.create(global_resolution, showWarnings = FALSE)

  xenium.obj@meta.data$cellname <- rownames(xenium.obj@meta.data)
  xenium.obj@meta.data = joindf(xenium.obj@meta.data,
    as.data.frame(xenium.obj@reductions$umap@cell.embeddings))

  filename_init <- object_file(opt$prefix)
  saveRDS(object = xenium.obj, file=filename_init)
  saveRDS(object = xenium.obj@reductions, file=gsub("init", "reductions", filename_init))
  saveRDS(object = xenium.obj@meta.data, file=gsub("init", "metadata", filename_init))
  saveRDS(object = xenium.obj@graphs, file=gsub("init", "graphs", filename_init))
}

if(opt$stage == "normalize_base") {
  run_normalize_base()
} else if(opt$stage == "hvg") {
  run_hvg()
} else if(opt$stage == "pca") {
  run_pca()
} else if(opt$stage == "cluster") {
  run_cluster()
} else {
  stop("Unknown --stage: ", opt$stage)
}

if(opt$verbose){
  cat("\n\n*******************************************************************\n")
  cat("Starting time:\n"); cat(st.time, "\n")
  cat("Finishing time:\n"); timestamp()
  cat("*******************************************************************\n")
  cat("SESSION INFO:\n"); print(sessionInfo()); cat("\n")
  cat("Pipeline finished successfully\n")
}
