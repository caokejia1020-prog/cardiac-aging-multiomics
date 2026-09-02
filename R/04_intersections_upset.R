formal_feature_set <- function(comparison_result) {
  unique(comparison_result$formal$Feature_ID)
}

build_intersection_for_sex <- function(layer_results, sex, cfg) {
  comparison_rows <- cfg$comparisons[cfg$comparisons$sex == sex, , drop = FALSE]
  comparison_ids <- comparison_rows$comparison_id
  sets <- lapply(comparison_ids, function(id) formal_feature_set(layer_results[[id]]))
  names(sets) <- comparison_ids
  universe <- sort(unique(unlist(sets, use.names = FALSE)))
  membership <- data.frame(Feature_ID = universe, stringsAsFactors = FALSE)
  for (id in comparison_ids) membership[[id]] <- universe %in% sets[[id]]

  if (!length(universe)) {
    counts <- data.frame(
      Pattern = "000",
      Intersection = "No formal differential features",
      Count = 0L,
      stringsAsFactors = FALSE
    )
  } else {
    bool <- as.matrix(membership[, comparison_ids, drop = FALSE])
    pattern <- apply(bool, 1, function(z) paste(as.integer(z), collapse = ""))
    membership$Pattern <- pattern
    count_table <- sort(table(pattern), decreasing = TRUE)
    counts <- data.frame(Pattern = names(count_table), Count = as.integer(count_table),
                         stringsAsFactors = FALSE)
    short_labels <- paste0(comparison_rows$first_label, "-", comparison_rows$second_label)
    counts$Intersection <- vapply(counts$Pattern, function(code) {
      bits <- strsplit(code, "", fixed = TRUE)[[1]] == "1"
      paste(short_labels[bits], collapse = " & ")
    }, character(1))
  }
  set_sizes <- data.frame(
    Sex = sex,
    Comparison = comparison_ids,
    Set_size = vapply(sets, length, integer(1)),
    stringsAsFactors = FALSE
  )
  counts$Sex <- sex
  list(membership = membership, counts = counts, set_sizes = set_sizes,
       comparison_rows = comparison_rows)
}

build_all_intersections <- function(all_results, cfg) {
  out <- list()
  for (layer_id in names(all_results)) {
    out[[layer_id]] <- list(
      Female = build_intersection_for_sex(all_results[[layer_id]], "Female", cfg),
      Male = build_intersection_for_sex(all_results[[layer_id]], "Male", cfg)
    )
  }
  out
}

plot_upset_sex <- function(intersection, sex) {
  counts <- intersection$counts
  comparison_rows <- intersection$comparison_rows
  short_labels <- paste0(comparison_rows$first_label, "-", comparison_rows$second_label)
  comparison_ids <- comparison_rows$comparison_id
  counts <- counts[order(-counts$Count, counts$Pattern), , drop = FALSE]
  counts$Intersection <- factor(counts$Intersection, levels = counts$Intersection)

  bar_plot <- ggplot2::ggplot(counts, ggplot2::aes(x = .data$Intersection, y = .data$Count)) +
    ggplot2::geom_col(fill = if (sex == "Female") "#B34D8C" else "#3C78B5", width = 0.72) +
    ggplot2::geom_text(ggplot2::aes(label = .data$Count), vjust = -0.25, size = 3.1) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.14))) +
    ggplot2::labs(title = sex, x = NULL, y = "Intersection size") +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0)
    )

  matrix_rows <- list()
  for (i in seq_len(nrow(counts))) {
    code <- strsplit(as.character(counts$Pattern[i]), "", fixed = TRUE)[[1]]
    if (length(code) != length(comparison_ids)) code <- rep("0", length(comparison_ids))
    matrix_rows[[i]] <- data.frame(
      Intersection = counts$Intersection[i],
      Comparison = factor(short_labels, levels = rev(short_labels)),
      Present = code == "1",
      stringsAsFactors = FALSE
    )
  }
  matrix_data <- dplyr::bind_rows(matrix_rows)
  matrix_plot <- ggplot2::ggplot(
    matrix_data,
    ggplot2::aes(x = .data$Intersection, y = .data$Comparison)
  ) +
    ggplot2::geom_point(color = "#D9D9D9", size = 3.0) +
    ggplot2::geom_point(data = matrix_data[matrix_data$Present, , drop = FALSE],
                        color = "#222222", size = 3.2) +
    ggplot2::geom_line(
      data = matrix_data[matrix_data$Present, , drop = FALSE],
      ggplot2::aes(group = .data$Intersection), color = "#222222", linewidth = 0.7
    ) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 35, hjust = 1, vjust = 1),
      axis.text.y = ggplot2::element_text(color = "black")
    )
  patchwork::wrap_plots(bar_plot, matrix_plot, ncol = 1, heights = c(2.2, 1.0))
}

save_plot_formats <- function(plot, stem, cfg, width = NULL, height = NULL) {
  width <- width %||% cfg$figures$width
  height <- height %||% cfg$figures$height
  ensure_dir(dirname(stem))
  files <- character()
  for (format in cfg$figures$formats) {
    path <- paste0(stem, ".", format)
    if (format == "svg") {
      ggplot2::ggsave(path, plot, device = svglite::svglite, width = width, height = height,
                      units = "in", bg = "white")
    } else if (format == "png") {
      ggplot2::ggsave(path, plot, width = width, height = height, units = "in",
                      dpi = cfg$figures$dpi, bg = "white")
    } else {
      ggplot2::ggsave(path, plot, width = width, height = height, units = "in",
                      bg = "white")
    }
    files <- c(files, path)
  }
  files
}

plot_all_upsets <- function(intersections, cfg, start_panel = 25L) {
  figure_dir <- file.path(cfg$paths$output_dir, "figures", "main_28_panels")
  ensure_dir(figure_dir)
  manifest <- list()
  for (i in seq_along(cfg$layers$layer_id)) {
    layer_id <- cfg$layers$layer_id[i]
    display <- cfg$layers$display_name[i]
    female <- plot_upset_sex(intersections[[layer_id]]$Female, "Female")
    male <- plot_upset_sex(intersections[[layer_id]]$Male, "Male")
    combined <- patchwork::wrap_plots(female, male, ncol = 1) + patchwork::plot_annotation(
      title = paste0(display, ": shared and unique formal limma differences"),
      subtitle = "Formal rule: BH-FDR < 0.05 and |log2FC| >= log2(1.5)",
      theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 14))
    )
    panel_no <- start_panel + i - 1L
    stem <- file.path(figure_dir, sprintf("P%02d_%s_SharedUnique_UpSet", panel_no, layer_id))
    files <- save_plot_formats(combined, stem, cfg, width = 9.0, height = 9.0)
    manifest[[length(manifest) + 1L]] <- data.frame(
      Panel = sprintf("P%02d", panel_no),
      Panel_type = "Shared_unique_UpSet",
      Layer = layer_id,
      Comparison = "Female and Male age contrasts",
      Primary_file = files[grepl("\\.png$", files)],
      stringsAsFactors = FALSE
    )
  }
  dplyr::bind_rows(manifest)
}
