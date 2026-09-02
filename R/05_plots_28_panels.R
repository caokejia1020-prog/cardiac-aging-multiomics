feature_label_for_plot <- function(result, label_column) {
  label <- if (label_column %in% names(result)) result[[label_column]] else result$Feature_ID
  label <- trim_na(label)
  label <- sub(";.*$", "", label)
  label[is.na(label)] <- result$Feature_ID[is.na(label)]
  label
}

plot_volcano <- function(result, layer_display, comparison_row, label_column, cfg) {
  plot_data <- result[result$Testable & is.finite(result$logFC) & is.finite(result$adj.P.Val), , drop = FALSE]
  plot_data$Plot_call <- factor(as.character(plot_data$Formal_call),
                               levels = c("Up", "Down", "Not_significant"))
  plot_data$minus_log10_FDR <- safe_neg_log10(plot_data$adj.P.Val)
  plot_data$Label <- feature_label_for_plot(plot_data, label_column)
  formal <- plot_data[plot_data$Plot_call %in% c("Up", "Down"), , drop = FALSE]
  label_rows <- formal |>
    dplyr::group_by(.data$Plot_call) |>
    dplyr::arrange(.data$adj.P.Val, dplyr::desc(abs(.data$logFC)), .by_group = TRUE) |>
    dplyr::slice_head(n = cfg$figures$volcano_label_n_each_direction) |>
    dplyr::ungroup()

  x_limit <- max(abs(plot_data$logFC), cfg$analysis$formal$minimum_abs_log2FC, na.rm = TRUE)
  x_limit <- ceiling(x_limit * 1.08 * 2) / 2
  title <- sprintf("%s: %s %s vs %s", layer_display, comparison_row$sex,
                   comparison_row$first_label, comparison_row$second_label)
  subtitle <- sprintf("Two-group limma; formal: BH-FDR < %.2g and |log2FC| >= %.3f",
                      cfg$analysis$formal$p_cutoff,
                      cfg$analysis$formal$minimum_abs_log2FC)

  ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$logFC, y = .data$minus_log10_FDR)) +
    ggplot2::geom_point(ggplot2::aes(color = .data$Plot_call), alpha = 0.72, size = 1.45) +
    ggplot2::geom_vline(xintercept = c(-1, 1) * cfg$analysis$formal$minimum_abs_log2FC,
                        linetype = 2, color = "#4D4D4D", linewidth = 0.45) +
    ggplot2::geom_hline(yintercept = -log10(cfg$analysis$formal$p_cutoff),
                        linetype = 2, color = "#4D4D4D", linewidth = 0.45) +
    ggrepel::geom_text_repel(
      data = label_rows,
      ggplot2::aes(label = .data$Label),
      size = 2.7, min.segment.length = 0, max.overlaps = Inf,
      box.padding = 0.30, point.padding = 0.15, seed = cfg$project$seed,
      show.legend = FALSE
    ) +
    ggplot2::scale_color_manual(values = cfg$figures$colors, drop = FALSE) +
    ggplot2::coord_cartesian(xlim = c(-x_limit, x_limit)) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = sprintf("log2 fold change (%s - %s)", comparison_row$first_label,
                  comparison_row$second_label),
      y = expression(-log[10](BH-FDR)),
      color = NULL
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(linewidth = 0.25, color = "#E8E8E8"),
      plot.title = ggplot2::element_text(face = "bold"),
      legend.position = "top"
    )
}

plot_all_differential_panels <- function(all_results, cfg) {
  figure_dir <- file.path(cfg$paths$output_dir, "figures", "main_28_panels")
  ensure_dir(figure_dir)
  manifest <- list()
  panel_no <- 0L
  for (layer_index in seq_len(nrow(cfg$layers))) {
    layer_id <- cfg$layers$layer_id[layer_index]
    layer_display <- cfg$layers$display_name[layer_index]
    label_column <- cfg$layers$label_column[layer_index]
    for (comparison_index in seq_len(nrow(cfg$comparisons))) {
      panel_no <- panel_no + 1L
      comparison <- cfg$comparisons[comparison_index, , drop = FALSE]
      result <- all_results[[layer_id]][[comparison$comparison_id]]$all
      plot <- plot_volcano(result, layer_display, comparison, label_column, cfg)
      stem <- file.path(
        figure_dir,
        sprintf("P%02d_%s_%s", panel_no, layer_id, comparison$comparison_id)
      )
      files <- save_plot_formats(plot, stem, cfg)
      manifest[[length(manifest) + 1L]] <- data.frame(
        Panel = sprintf("P%02d", panel_no),
        Panel_type = "Differential_volcano",
        Layer = layer_id,
        Comparison = comparison$comparison_id,
        Primary_file = files[grepl("\\.png$", files)],
        stringsAsFactors = FALSE
      )
    }
  }
  if (panel_no != cfg$expected$number_of_differential_panels) {
    abortf("Expected %d differential panels, generated %d",
           cfg$expected$number_of_differential_panels, panel_no)
  }
  dplyr::bind_rows(manifest)
}
