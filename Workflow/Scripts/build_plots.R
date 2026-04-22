suppressPackageStartupMessages({
  library(tidyverse)
  library(ggpubr)
  library(patchwork)
  library(scales)
})


#----- Parameters from snakemake@params -----#
plot_type     <- snakemake@params$type
basesize      <- as.integer(snakemake@params$basesize %||% 12)
pvalue_cutoff <- as.numeric(snakemake@params$pvalue   %||% 0.01)
locus_tags    <- snakemake@params$locus %||% character(0)
gene_names    <- snakemake@params$gene  %||% character(0)


#----- Read each variant CSV, tag with a Calibration column based on its path -----#
read_one <- function(path) {
  variant <- if (grepl("/w_calibration/", path, fixed = TRUE)) "calibrated" else "uncalibrated"
  read_csv(path, show_col_types = FALSE) |> mutate(Calibration = variant)
}

merged <- map(snakemake@input$annotated, read_one) |>
  bind_rows() |>
  mutate(Calibration = factor(Calibration, levels = c("uncalibrated", "calibrated")))


#----- Validated-hits filter (locus_tag OR gene_name match), one row per (Calibration, miRNA, locus_tag) -----#
validated <- merged |>
  filter(locus_tag %in% locus_tags | gene_name %in% gene_names) |>
  group_by(Calibration) |>
  arrange(P_value) |>
  distinct(miRNA, locus_tag, .keep_all = TRUE) |>
  ungroup()


#----- Auto color + label maps for the validated hits, keyed by locus_tag -----#
hit_levels <- unique(validated$locus_tag)
hit_colors <- setNames(scales::hue_pal()(length(hit_levels)), hit_levels)
hit_label_table <- validated |>
  distinct(locus_tag, miRNA, gene_name) |>
  mutate(label = sprintf("%s | %s", miRNA, coalesce(gene_name, locus_tag)))
hit_labels <- setNames(hit_label_table$label, hit_label_table$locus_tag)
hit_shapes <- setNames(seq_along(hit_levels) + 14L, hit_levels)


#----- Pretty labels for the Calibration facet -----#
facet_labels <- as_labeller(c(
  "uncalibrated" = "Without Calibration",
  "calibrated"   = "With Calibration"
))


#----- Plot builders -----#
plot_pvalue_distribution <- function() {
  # Auto Y-axis cut: ceiling(max per-facet density y-max * 1.05)
  dens <- merged |>
    group_by(Calibration) |>
    summarise(ymax = max(density(P_value, na.rm = TRUE)$y), .groups = "drop")
  ycut <- ceiling(max(dens$ymax) * 1.05)

  ggplot(merged, aes(x = P_value)) +
    geom_density(color = "black", fill = "#d1d1d1") +
    geom_vline(
      data = validated,
      aes(xintercept = P_value, color = locus_tag),
      linetype = "dashed"
    ) +
    scale_color_manual(name = "Validated hits", values = hit_colors, labels = hit_labels) +
    coord_cartesian(ylim = c(0, ycut)) +
    facet_wrap(~Calibration, labeller = facet_labels) +
    labs(
      x = "p-value", y = "Density",
      title = sprintf("Distribution of p-values (Y-axis %d)", ycut)
    ) +
    theme_bw(base_size = basesize) +
    theme(strip.placement = "outside")
}

plot_position_pvalue <- function() {
  ggplot(merged, aes(x = Position, y = -log10(P_value))) +
    geom_point(
      data = merged |> filter(P_value <= pvalue_cutoff),
      color = "black", size = 0.7
    ) +
    geom_point(
      data = validated,
      aes(shape = locus_tag),
      color = "red"
    ) +
    scale_shape_manual(name = "Validated hits", values = hit_shapes, labels = hit_labels) +
    geom_smooth(method = "lm", linewidth = 0.5) +
    stat_cor(method = "spearman") +
    labs(
      x = "Relative Position", y = "-log10(p-value)",
      title = sprintf("Position vs p-values (≤ %g)", pvalue_cutoff)
    ) +
    theme_bw(base_size = basesize) +
    facet_wrap(~Calibration, labeller = facet_labels) +
    theme(strip.placement = "outside")
}

plot_position_energy <- function() {
  ggplot(merged, aes(x = Position, y = Energy)) +
    geom_point(
      data = merged |> filter(P_value <= pvalue_cutoff),
      color = "black", size = 0.7
    ) +
    geom_point(
      data = validated,
      aes(shape = locus_tag),
      color = "red"
    ) +
    scale_shape_manual(name = "Validated hits", values = hit_shapes, labels = hit_labels) +
    geom_smooth(method = "lm", linewidth = 0.5) +
    stat_cor(method = "spearman") +
    labs(
      x = "Relative Position", y = "Energy (KCal/Mol)",
      title = sprintf("Position vs Energy (p-value ≤ %g)", pvalue_cutoff)
    ) +
    facet_wrap(~Calibration, labeller = facet_labels) +
    theme_bw(base_size = basesize) +
    theme(strip.placement = "outside")
}


#----- Dispatch on type -----#
plots <- switch(plot_type,
  "pvalue_distribution" = list(plot_pvalue_distribution()),
  "position_pvalue"     = list(plot_position_pvalue()),
  "position_energy"     = list(plot_position_energy()),
  "all"                 = list(
    plot_pvalue_distribution(),
    plot_position_pvalue(),
    plot_position_energy()
  ),
  stop(sprintf("Unknown build_plots.type: %s", plot_type))
)


#----- Compose & save -----#
combined <- if (length(plots) == 1) {
  plots[[1]]
} else {
  Reduce(`|`, plots) + plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
}

dir.create(dirname(snakemake@output$pdf), showWarnings = FALSE, recursive = TRUE)
ggsave(
  snakemake@output$pdf, combined,
  width = 6 * length(plots), height = 5, units = "in"
)
