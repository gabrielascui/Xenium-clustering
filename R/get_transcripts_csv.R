#!/usr/bin/R

# ---
# Author: Francisco Emmanuel Castaneda-Castro
# Other contributors: 
    # Donaldo Sosa-Garcia
# Taking as a base Ciro's single-cell clustering pipeline (https://github.com/vijaybioinfo/clustering)
# Date: 2025-06-24
# Version 1.0.0
## This script will only process cellranger output and create the transcript.csv file needed. 
    # v2.0.0: Customization for multiple slides

# ---

######################
# Clustering: Seurat: Spatial analysis # Step 1: get_transcripts_csv
######################
lib4.3path = "/home/fcastaneda/R/x86_64-pc-linux-gnu-library/4.3"
if(grepl("4.3", getRversion()) && file.exists(lib4.3path))
  .libPaths(new = lib4.3path)
options(future.globals.maxSize= 13312*1024^2)

library(optparse)
optlist <- list(
  optparse::make_option(
    opt_str = c("-y", "--yaml"), type = "character",
    help = "Configuration file: Instructions in YAML format."
  ),
    optparse::make_option(
    opt_str = c("-v", "--verbose"), default = TRUE,
    help = "Verbose: Show progress."
  ),
    optparse::make_option(
    opt_str = c("-s", "--slide"), default = 1,
    help = "Verbose: Show progress."
  )
)

# Getting arguments from command line
opt <- optparse::parse_args(optparse::OptionParser(option_list = optlist))

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if(length(script_file) && !is.na(script_file)) dirname(normalizePath(script_file)) else "R"
source(file.path(script_dir, "functions.R"))

get_transcripts_csv<-function(cellranger_path){ 
    cat ("transcripts.csv file not present \n"); 
    cat ("Converting transcripts.parquet to csv as required for Seurat \n"); 
    # cat ("Creating transcripts files \n"); 
      ## Thanks cellranger https://www.10xgenomics.com/support/software/xenium-onboard-analysis/latest/advanced/example-code 13 Nov 2024
      # cellranger_path <- cellranger_out
        library(arrow)
        # Path to your parquet file, edit path to where parquet file saved
        PATH <- paste0(cellranger_path, '/transcripts.parquet')
        # Edit path and output name for new file
        OUTPUT <- gsub('\\.parquet$', '.csv', PATH)
        # Specify chunk size
        CHUNK_SIZE <- 1e6

        cat("Read in the parquet file \n")
        # This reads in the data as an arrow table
        parquet_file <- arrow::read_parquet(PATH, as_data_frame = FALSE)
        # To read in the table as a tibble, set data frame to true:
        # parquet_file <- arrow::read_parquet(PATH, as_data_frame = TRUE)
        #Optional: convert parquet data frame to CSV.

        cat("Writting: convert parquet data frame to CSV. \n") 
        start <- 0
        while(start < parquet_file$num_rows) {
          end <- min(start + CHUNK_SIZE, parquet_file$num_rows)
          chunk <- as.data.frame(parquet_file$Slice(start, end - start))
          data.table::fwrite(chunk, OUTPUT, append = start != 0)
          start <- end
        }
        if(require('R.utils', quietly = TRUE)) {
          R.utils::gzip(OUTPUT)
        }

}

# opt parameters have the priority
if(interactive()){ # Example/manually
  opt$yaml = "config.yaml"
  opt$verbose <- TRUE
  opt$slide <- 1
}

config = yaml::read_yaml(opt$yaml)
outdir_sp = "." #This need to change to an option to run it without the snakemake
if(interactive()) outdir_sp= paste0(config$output_dir, "/", config$project_name)

config = yaml::read_yaml(opt$yaml)
cellranger_out= config$input_expression[opt$slide]
setwd(outdir_sp)
cat("Working in:", getwd(), "\n")

if(opt$verbose) cat('Date and time:\n') ; st.time <- timestamp();
if(opt$verbose) cat("Current directory: ", getwd() , "\n") ;

if(grepl("init|.rds", cellranger_out)) { 
  cat("Assumming this is not the first clustering run \n ") 
  file_type <- "RDS_file_provided"
} else {
  file_type <- "parquet"
## Checking if parquet files are already in csv as required for Seurat
  if(!dir.exists(cellranger_out)) stop("Xenium-ranger result folder doesn't exist") # probably this is unneessary in the snakemake
  files_cellranger <- list.files(cellranger_out)
  if(!any(grepl("transcripts.csv.gz", files_cellranger))) { get_transcripts_csv(cellranger_out) } else cat("transcripts.csv.gz exists for slide ", opt$slide)
}
writeLines(paste0("transcripts.csv.gz from ",file_type ," file processed"), paste0("transcripts_done_slide", opt$slide, ".txt"))


if(opt$verbose){
  cat('\n\n*******************************************************************\n')
  cat('Starting time:\n'); cat(st.time, '\n')
  cat('Finishing time:\n'); timestamp()
  cat('*******************************************************************\n')
  cat('SESSION INFO:\n'); print(sessionInfo()); cat("\n")
  cat('Pipeline finished successfully\n')
}
