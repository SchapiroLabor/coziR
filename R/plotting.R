#' Dotplot for COZI Results
#'
#' Creates a dotplot where color encodes z-score and dot size encodes
#' conditional cell ratio (`con_cell_ratio`).
#'
#' @param results Output list from [run_cozi()] or [nep_analysis()], containing
#'   `zscore` and `cond_ratio`.
#' @param title Plot title.
#' @param low,mid,high Colors used for low/mid/high z-score values.
#' @param size_limits Numeric length-2 vector for conditional ratio size scale.
#'   Defaults to `c(0, 1)`.
#'
#' @return A `ggplot2` object.
#' @export
cozi_dotplot <- function(results,
                         title = "COZI Dotplot",
                         low = "#2166AC",
                         mid = "white",
                         high = "#B2182B",
                         size_limits = c(0, 1)) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package `ggplot2` is required for `cozi_dotplot()`.", call. = FALSE)
  }

  if (!is.list(results) || !all(c("zscore", "cond_ratio") %in% names(results))) {
    stop("`results` must be a list containing `zscore` and `cond_ratio`.", call. = FALSE)
  }

  zscore <- as.matrix(results$zscore)
  cond_ratio <- as.matrix(results$cond_ratio)

  if (!all(dim(zscore) == dim(cond_ratio))) {
    stop("`zscore` and `cond_ratio` must have matching dimensions.", call. = FALSE)
  }

  row_ids <- rownames(zscore)
  col_ids <- colnames(zscore)
  if (is.null(row_ids)) {
    row_ids <- as.character(seq_len(nrow(zscore)))
  }
  if (is.null(col_ids)) {
    col_ids <- as.character(seq_len(ncol(zscore)))
  }

  plot_df <- data.frame(
    source_cell_type = rep.int(row_ids, times = ncol(zscore)),
    neighbor_cell_type = rep(col_ids, each = nrow(zscore)),
    zscore = as.numeric(zscore),
    con_cell_ratio = as.numeric(cond_ratio),
    stringsAsFactors = FALSE
  )

  ggplot2::ggplot(
    plot_df,
    ggplot2::aes_string(
      x = "neighbor_cell_type",
      y = "source_cell_type",
      color = "zscore",
      size = "con_cell_ratio"
    )
  ) +
    ggplot2::geom_point(alpha = 0.9) +
    ggplot2::scale_color_gradient2(
      low = low,
      mid = mid,
      high = high,
      midpoint = 0,
      name = "z-score"
    ) +
    ggplot2::scale_size_continuous(
      name = "con_cell_ratio",
      limits = size_limits
    ) +
    ggplot2::labs(
      x = "Neighbor cell type",
      y = "Source cell type",
      title = title
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
}
