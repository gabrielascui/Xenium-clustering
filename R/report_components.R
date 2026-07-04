#!/usr/bin/R

# ---
# Author: Francisco Emmanuel Castaneda-Castro
# Other contributors: 
    # Donaldo Sosa-Garcia
# Taking as a base Ciro's single-cell clustering pipeline (https://github.com/vijaybioinfo/clustering)
# Date: 2025-06-04
# Version 1.0.0
    # This code does the visualization of the dimensionality reduction. 
# ---


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
    opt_str = c("-n", "--prefix"), type = "character",
    help = "Prefix for output. Name of the files."
  ),
  optparse::make_option(
    opt_str = c("-i", "--init_file"), type = "character",
    help = "Prefix for output. Name of the files."
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

ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

# opt parameters have the priority
if(interactive()){ # Example/manually
  opt$yaml = "config.yaml"
  opt$verbose <- TRUE
  opt$percent <- 30
  opt$chosen_comp <- 15
  opt$init_file <- ".object_init_mean0.01_pct30_pc15.rds"
  opt$prefix <- "init_mean0.01_pct30_pc15"
}

config = yaml::read_yaml(opt$yaml)
cellranger_out= config$input_expression
pcts= opt$percent
pc= opt$chosen_comp
if(interactive()) outdir_sp= paste0(config$output_dir, "/", config$project_name)
outdir_sp = "."

if(opt$verbose) cat('Date and time:\n') ; st.time <- timestamp();
## Chekinf if parquet files are already in csv as required for Seurat
if(!dir.exists(outdir_sp)) dir.create(outdir_sp, recursive=TRUE)
setwd(outdir_sp)
cat("Working in:", getwd(), "\n")
# Load the init object --> Xenium

checking_init <- paste0("pct", pcts, "_pc", pc, ".rds")
if(!grepl(checking_init, opt$init_file)) stop("Init file doesn't mach with pct or pc given!")
xenium.obj <- readRDS(opt$init_file)
resolutions <- config$resolution
  lapply(resolutions, function(resolution_cho){
        # resolution_cho<-resolutions[1]
        comn_resolution<<-paste0(pcts, "pct", pc, "pc_res",resolution_cho)
        if(!dir.exists(comn_resolution)) dir.create(comn_resolution)
        cat(paste("Performing: ",comn_resolution, "\n"))
        cluster_value<- ifelse(grepl(config$norm, "sctransform"), paste0("SCT_snn_res.", resolution_cho), paste0("Xenium_snn_res.", resolution_cho))

        xenium.obj@meta.data[, cluster_value] <- factor(xenium.obj@meta.data[, cluster_value], levels=c(0:length(unique(xenium.obj@meta.data[, cluster_value]))))

        # VariableFeatures(object = xenium.obj) 

        pdf(paste0(comn_resolution,"/1.dimentional_reduction_umap_res", resolution_cho, ".pdf"))
        print(DimPlot(xenium.obj, group.by=cluster_value, label =T ))
        dev.off()

        umap_plot <- DimPlot(xenium.obj, group.by=cluster_value, size = 2, label =T)
        interactive_plot <- ggplotly(umap_plot)

        # Save as HTML
        htmlwidgets::saveWidget(interactive_plot, paste0(comn_resolution,"/1.dimentional_reduction_umap_res", resolution_cho, ".html"))

        pdf(paste0(comn_resolution, "/2.qc_violins_per_cluster.pdf"), width=8)
        aver<-VlnPlot(xenium.obj, features = c("nFeature_Xenium", "nCount_Xenium"), group.by=cluster_value,ncol = 2, pt.size=0.00, alpha=0.2, combine=FALSE)
        print(aver[[1]] + geom_boxplot(width=0.2))
        print(aver[[2]] + geom_boxplot(width=0.2))
        dev.off()

        pdf(paste0(comn_resolution, "/3.UMAP_inSpatial.pdf"), height=36, width=84)
        for (fov in names(xenium.obj@images)){
        print(ImageDimPlot(xenium.obj, size = 2, fov = fov, border.size=NA,dark.background = TRUE, group.by=cluster_value))
        }
        dev.off() 

        umap_tissue_plot <- ImageDimPlot(xenium.obj, size = 0.5, border.size=NA,dark.background = TRUE, group.by=cluster_value)

        interactive_plot <- ggplotly(umap_tissue_plot)
        htmlwidgets::saveWidget(interactive_plot, paste0(comn_resolution,"/3.UMAP_inSpatial", ".html"))

        pdf(paste0(comn_resolution, "/3.UMAP_inSpatial_per_cluster.pdf"), height=54, width=126)
        for (fov in names(xenium.obj@images)){
        print(ImageDimPlot(xenium.obj, size = 5, fov = fov, border.size=NA,dark.background = TRUE, group.by=cluster_value, split.by = cluster_value))
        }
        dev.off() 

        xenium.obj@meta.data$cellnames <- rownames(xenium.obj@meta.data)
        cells<-data.frame(
          cell_id =xenium.obj$cellnames,
          group = xenium.obj@meta.data[, cluster_value]
        )

        write.csv(cells, paste0(comn_resolution, "/5.to_xenium_explorer", comn_resolution, "_cell_class.csv"), row.names=FALSE)
        
      cat("Plotting proportions \n")

        colrs<-sort(as.numeric(unique(xenium.obj@meta.data[, cluster_value])))
        lcolrs<-length(colrs) + 1
        color_vector<-ggplotColours(n=lcolrs)[-1]
        names(color_vector) <- colrs
        # cat("Using random colors: \n", names(color_vector), "\n")

          gcols <- color_vector
          confounders<- grep("orig.", colnames(xenium.obj@meta.data), value=TRUE)
          metadata <- xenium.obj@meta.data
          resdir <- paste0(comn_resolution, "/")
          prop.normalise = TRUE

          confounders <- grep("V1|cells", confounders, value=TRUE, invert=TRUE)

          metadata$Identity <- factor(
          metadata[, cluster_value],
          levels = mixedsort(unique(as.character(metadata[, cluster_value]))))


        centers <- list()
        reductions <- ident_dimensions(colnames(metadata), "umap_")
        for(rdim in names(reductions)){
          metadata$x1234 <- metadata[, reductions[[rdim]][1]];
          metadata$y1234 <- metadata[, reductions[[rdim]][2]]
          centers[[rdim]] <- metadata %>% group_by(Identity) %>%
            summarize(x = median(x = x1234), y = median(x = y1234))
        }; metadata <- metadata[, !colnames(metadata) %in% c("x1234", "y1234")]

        if(length(levels(metadata$Identity)) > 1) {
          for(orig in confounders){
            # orig<-confounders[1]
            if(opt$verbose) cat("   #", orig, "\n")
            df_plots <- metadata[!is.na(metadata[, orig]), ]
            df_plots <- df_plots[df_plots[, orig] != "", ]
            identity_levels <- mixedsort(unique(as.character(df_plots[, orig])))
            df_plots$var_interest <- factor(df_plots[, orig], levels = identity_levels)
            identities <- table(df_plots$var_interest) # needs the size per identity
            if(length(identities) == 1){ if(opt$verbose) cat(" (N = 1)\n"); next }
            fname <- paste0(resdir, 'proportions_', sub("orig\\.|orig", "", orig))
            fnames <- paste0(
              fname, c(".csv", paste0("_", names(reductions), ".pdf"))
            )
            if(all(is.file.finished(fnames))){ if(opt$verbose) cat(" done\n"); next }
            if(opt$verbose) cat("\n")

            do_i_norm <- !grepl("ht_id", orig)
            tvar <- (min(identities) < 500) && do_i_norm
            thesecells <- if(((min(identities) / sum(identities)) < 0.1) && tvar){
              if(opt$verbose)
                cat("     Minimum size cluster is < 10 %, sampling to a median of",
                  median(identities), "per group\n")
              sample_even(annot = df_plots, cname = orig,
                maxln = -median(identities), v = FALSE)
            }else if(do_i_norm){
              if(opt$verbose) cat("     Random sampling\n")
              sample_even(annot = df_plots, cname = orig, v = TRUE, maxln=20000)
            }else{ rownames(df_plots) }

            thesecells <- unname(thesecells)
            p_bar <- plot_pct(
              x = df_plots, groups = c("var_interest", "Identity"),
              normalise = prop.normalise && do_i_norm, return_table = TRUE, v=TRUE
            )
            write.csv(p_bar$table, file = fnames[1], quote = FALSE)
            p_bar <- p_bar$plot +
              guides(fill = list(ncol = make_grid(nlevels(df_plots$var_interest))[2])) +
              scale_fill_manual(values = v2cols(names(identities), gcols)) + coord_flip();
            p_pie <- plot_pct(
              x = df_plots, groups = c("Identity", "var_interest"),
              normalise = FALSE, type = "pie") +
              scale_fill_manual(values = v2cols(levels(df_plots$Identity), gcols)) +
              labs(fill = NULL) + guides(colour = guide_legend(override.aes = list(size = 6)))

            if(opt$verbose) cat("     dimensional reductions \n")
            print("hola /n")
            print(names(reductions))
            for(rdim in names(reductions)){ # Visualising with sampled data. For loop is not neccesary. Keeping just in case we anated to add another dimensionality reduction
              # rdim <- names(reductions)[1]
              fname_i = paste0(fname, "_", rdim, ".pdf"); if(opt$verbose) cat(" -", rdim, "\n")
              if(is.file.finished(fname_i)) next
              pp <- plot_grids(df_plots[thesecells, ], x = reductions[[rdim]][1],
                y = reductions[[rdim]][2], color = "Identity", colours = gcols,
                centers = centers[[rdim]],
                facet = list("var_interest", make_grid(identities)[1]))
              pp2 <- plot_grids(df_plots[thesecells, ], x = reductions[[rdim]][1],
                y = reductions[[rdim]][2], color = "var_interest",
                centers = centers[[rdim]])
              if(nlevels(df_plots$var_interest) == 2)
                pp <- pp + theme(plot.margin = margin(0, 3, 0, 3, "cm"))
              pdf(fname_i, width = 15, height = 10)
              print(cowplot::plot_grid(pp, p_bar + NoLegend(), rel_widths = c(2, 0.8), ncol = 2))
              print(cowplot::plot_grid(pp2, p_pie + NoLegend(), rel_widths = c(1.5, 0.8), ncol = 2))
              graphics.off()

              fname_i = paste0(fname, "_", rdim, "_inSpatial.pdf"); if(opt$verbose) cat(" -", rdim, "\n")
              
              pdf(fname_i, height=54, width=126)
              for (fov in names(xenium.obj@images)){
              print(ImageDimPlot(xenium.obj, size = 5, fov = fov, border.size=NA,dark.background = TRUE, group.by=orig))
              }
              dev.off() 
            }; cat("\n")



          }
        }

            rm(xenium.obj);
  })

if(opt$verbose) cat("Saving output file for snakemake \n")
done_file <- paste0(".", gsub("init", "report", opt$prefix), ".txt")
writeLines("report components done", done_file)

if(opt$verbose){
  cat('\n\n*******************************************************************\n')
  cat('Starting time:\n'); cat(st.time, '\n')
  cat('Finishing time:\n'); timestamp()
  cat('*******************************************************************\n')
  cat('SESSION INFO:\n'); print(sessionInfo()); cat("\n")
  cat('Pipeline finished successfully\n')
}
  
