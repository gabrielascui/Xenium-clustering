## Local helper functions for the spatial clustering pipeline.
## This file replaces dependencies on lab-adjacent helper scripts under /home.

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

readfile <- function(file, ...) {
  if (!file.exists(file)) stop("File not found: ", file)
  if (grepl("\\.rds$", file, ignore.case = TRUE)) return(readRDS(file))
  if (grepl("\\.csv(\\.gz)?$", file, ignore.case = TRUE)) {
    return(read.csv(file, stringsAsFactors = FALSE, check.names = FALSE, ...))
  }
  if (grepl("\\.tsv(\\.gz)?$|\\.txt(\\.gz)?$", file, ignore.case = TRUE)) {
    return(read.delim(file, stringsAsFactors = FALSE, check.names = FALSE, ...))
  }
  read.table(file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, ...)
}

remove.factors <- function(x) {
  if (!is.data.frame(x)) return(x)
  x[] <- lapply(x, function(col) {
    if (is.factor(col)) as.character(col) else col
  })
  x
}

joindf <- function(x, y) {
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  y <- as.data.frame(y, stringsAsFactors = FALSE)
  common <- intersect(rownames(x), rownames(y))
  x <- x[common, , drop = FALSE]
  y <- y[common, , drop = FALSE]
  y <- y[, !colnames(y) %in% colnames(x), drop = FALSE]
  cbind(x, y)
}

is.file.finished <- function(files) {
  file.exists(files) & file.info(files)$size > 0
}

mixedsort <- function(x) {
  if (requireNamespace("gtools", quietly = TRUE)) return(gtools::mixedsort(x))
  sort(x)
}

make_grid <- function(n) {
  n <- length(n)
  if (n < 1) return(c(1, 1))
  nr <- ceiling(sqrt(n))
  nc <- ceiling(n / nr)
  c(nr, nc)
}

v2cols <- function(values, colours) {
  values <- as.character(values)
  if (is.null(names(colours))) {
    out <- rep(colours, length.out = length(values))
    names(out) <- values
    return(out)
  }
  missing_values <- setdiff(values, names(colours))
  if (length(missing_values) > 0) {
    extra <- grDevices::hcl.colors(length(missing_values), "Dark 3")
    names(extra) <- missing_values
    colours <- c(colours, extra)
  }
  colours[values]
}

filters_complex <- function(mdata, filters, verbose = TRUE) {
  keep <- rep(TRUE, nrow(mdata))
  for (flt in filters) {
    fname <- flt[[1]]
    expr <- flt[[2]]
    this_keep <- eval(parse(text = expr), envir = mdata, enclos = parent.frame())
    this_keep[is.na(this_keep)] <- FALSE
    if (verbose) cat("Filter", fname, "keeps", sum(this_keep), "of", length(this_keep), "rows\n")
    keep <- keep & this_keep
  }
  list(mdata[keep, , drop = FALSE])
}

qc_violin_spatial <- function(dat, yax, xax = "origlib", to_text = NULL) {
  p <- ggplot2::ggplot(dat, ggplot2::aes(x = .data[[xax]], y = .data[[yax]], fill = .data[[xax]])) +
    ggplot2::geom_violin(scale = "width", trim = TRUE, alpha = 0.7, na.rm = TRUE) +
    ggplot2::geom_boxplot(width = 0.12, outlier.size = 0.2, alpha = 0.8, na.rm = TRUE) +
    ggplot2::theme_classic() +
    ggplot2::labs(x = xax, y = yax, fill = xax)
  if (!is.null(to_text) && all(c(xax, "size") %in% colnames(to_text))) {
    p <- p + ggplot2::labs(subtitle = paste(paste(to_text[[xax]], to_text$size, sep = ": "), collapse = " | "))
  }
  p
}

VariableFeaturePlotp_bicolor <- function(object, assay = NULL, variable_features = NULL) {
  assay <- assay %||% Seurat::DefaultAssay(object)
  hvf <- Seurat::HVFInfo(object, assay = assay)
  hvf$gene <- rownames(hvf)
  variable_features <- variable_features %||% Seurat::VariableFeatures(object, assay = assay)
  hvf$variable <- hvf$gene %in% variable_features
  mean_col <- grep("^mean$|mean", colnames(hvf), value = TRUE)[1]
  disp_col <- grep("variance.standardized|dispersion.standardized|standardized", colnames(hvf), value = TRUE)[1]
  ggplot2::ggplot(hvf, ggplot2::aes(x = .data[[mean_col]], y = .data[[disp_col]], color = variable)) +
    ggplot2::geom_point(size = 0.7, alpha = 0.7) +
    ggplot2::scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "red3")) +
    ggplot2::theme_classic() +
    ggplot2::labs(x = mean_col, y = disp_col, color = "Variable")
}

sample_even <- function(annot, cname, maxln = 20000, v = TRUE) {
  groups <- split(rownames(annot), annot[[cname]])
  if (maxln < 0) maxln <- abs(maxln)
  sampled <- lapply(groups, function(cells) {
    n <- min(length(cells), maxln)
    if (length(cells) > n) sample(cells, n) else cells
  })
  unname(unlist(sampled))
}

ident_dimensions <- function(cols, pattern = "umap_") {
  hits <- grep(pattern, cols, value = TRUE, ignore.case = TRUE)
  if (length(hits) < 2) return(list())
  out <- list()
  for (i in seq(1, length(hits) - 1, by = 2)) {
    out[[gsub("[^A-Za-z0-9]+$", "", pattern) %||% "reduction"]] <- hits[i:(i + 1)]
  }
  out
}

plot_pct <- function(x, groups, normalise = TRUE, return_table = FALSE, type = "bar", v = TRUE) {
  df <- as.data.frame(x)
  tab <- as.data.frame(table(df[[groups[1]]], df[[groups[2]]]), stringsAsFactors = FALSE)
  colnames(tab) <- c(groups[1], groups[2], "n")
  if (normalise) {
    tab <- tab |>
      dplyr::group_by(.data[[groups[1]]]) |>
      dplyr::mutate(pct = ifelse(sum(n) > 0, n / sum(n) * 100, 0)) |>
      dplyr::ungroup()
    yvar <- "pct"
  } else {
    tab$pct <- tab$n
    yvar <- "n"
  }
  if (type == "pie") {
    p <- ggplot2::ggplot(tab, ggplot2::aes(x = "", y = .data[[yvar]], fill = .data[[groups[2]]])) +
      ggplot2::geom_col(width = 1) +
      ggplot2::coord_polar(theta = "y") +
      ggplot2::facet_wrap(stats::as.formula(paste("~", groups[1]))) +
      ggplot2::theme_void()
  } else {
    p <- ggplot2::ggplot(tab, ggplot2::aes(x = .data[[groups[1]]], y = .data[[yvar]], fill = .data[[groups[2]]])) +
      ggplot2::geom_col() +
      ggplot2::theme_classic() +
      ggplot2::labs(x = groups[1], y = ifelse(normalise, "Percent", "Count"))
  }
  if (return_table) return(list(plot = p, table = tab))
  p
}

plot_grids <- function(df, x, y, color, colours = NULL, centers = NULL, facet = NULL) {
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[x]], y = .data[[y]], color = .data[[color]])) +
    ggplot2::geom_point(size = 0.15, alpha = 0.7) +
    ggplot2::theme_classic() +
    ggplot2::labs(color = color)
  if (!is.null(colours)) p <- p + ggplot2::scale_color_manual(values = v2cols(unique(as.character(df[[color]])), colours))
  if (!is.null(centers) && all(c("x", "y") %in% colnames(centers))) {
    p <- p + ggplot2::geom_text(data = centers, ggplot2::aes(x = x, y = y, label = Identity),
      inherit.aes = FALSE, size = 3)
  }
  if (!is.null(facet)) p <- p + ggplot2::facet_wrap(stats::as.formula(paste("~", facet[[1]])), ncol = facet[[2]])
  p
}

make_list <- function(x, colname, grouping = TRUE) {
  split(rownames(x), x[[colname]])
}

.row_means <- function(mat) {
  if (inherits(mat, "IterableMatrix")) {
    return(BPCells::matrix_stats(mat, row_stats = "mean")$row_stats["mean", ])
  }
  Matrix::rowMeans(mat)
}

.row_nonzero_pct <- function(mat) {
  if (inherits(mat, "IterableMatrix")) {
    n <- ncol(mat)
    return(BPCells::matrix_stats(mat > 0, row_stats = "mean")$row_stats["mean", ] * 100)
  }
  Matrix::rowMeans(mat > 0) * 100
}

stats_summary_table <- function(mat, groups, rnames = rownames(mat), datatype = "SeuratNormalized", verbose = TRUE) {
  rnames <- intersect(rnames, rownames(mat))
  out <- lapply(names(groups), function(cluster) {
    cells <- intersect(groups[[cluster]], colnames(mat))
    if (length(cells) == 0) return(NULL)
    submat <- mat[rnames, cells, drop = FALSE]
    data.frame(
      gene = rnames,
      cluster = cluster,
      avg_log = as.numeric(.row_means(submat)),
      pct = as.numeric(.row_nonzero_pct(submat)),
      stringsAsFactors = FALSE
    )
  })
  dplyr::bind_rows(out)
}

markers_summary <- function(marktab, annot, datavis, cluster_column, datatype = "SeuratNormalized", verbose = TRUE) {
  groups <- make_list(annot, cluster_column, grouping = TRUE)
  genes <- unique(marktab$gene)
  stats <- stats_summary_table(datavis, groups = groups, rnames = genes, datatype = datatype, verbose = verbose)
  names(stats)[names(stats) == "pct"] <- "pct_cluster"
  bg <- lapply(names(groups), function(cluster) {
    other_cells <- setdiff(colnames(datavis), groups[[cluster]])
    if (length(other_cells) == 0) return(NULL)
    submat <- datavis[genes, other_cells, drop = FALSE]
    data.frame(
      gene = genes,
      cluster = cluster,
      avg_other = as.numeric(.row_means(submat)),
      pct_other = as.numeric(.row_nonzero_pct(submat)),
      stringsAsFactors = FALSE
    )
  })
  bg <- dplyr::bind_rows(bg)
  out <- dplyr::left_join(marktab, stats, by = c("gene", "cluster"))
  out <- dplyr::left_join(out, bg, by = c("gene", "cluster"))
  if (!"avg_log" %in% colnames(out)) {
    avg_col <- grep("^avg_log|avg_diff|avg_log2FC", colnames(out), value = TRUE)[1]
    out$avg_log <- out[[avg_col]]
  }
  out$Dpct <- out$pct_cluster - out$pct_other
  out
}

marker_report <- function(markers_df, file = "_summary_stats.csv", verbose = TRUE, return_plot = TRUE) {
  if (!return_plot) return(invisible(NULL))
  p1 <- ggplot2::ggplot(markers_df, ggplot2::aes(x = Dpct, y = avg_log, color = factor(cluster))) +
    ggplot2::geom_point(alpha = 0.7) +
    ggplot2::theme_classic() +
    ggplot2::labs(color = "Cluster")
  p2 <- ggplot2::ggplot(markers_df, ggplot2::aes(x = factor(cluster), y = Dpct, fill = factor(cluster))) +
    ggplot2::geom_boxplot(outlier.size = 0.2) +
    ggplot2::theme_classic() +
    ggplot2::labs(x = "Cluster", fill = "Cluster")
  list(p1, p2)
}
