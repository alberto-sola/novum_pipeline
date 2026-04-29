suppressPackageStartupMessages({
  library(tidyverse)
  library(ggpubr)
  library(patchwork)
  library(scales)
})


#----- Parameters from snakemake@params -----#
plot_type      <- snakemake@params$type
basesize       <- as.integer(snakemake@params$basesize %||% 12)
pvalue_cutoff  <- as.numeric(snakemake@params$pvalue   %||% 0.01)
locus_tags     <- snakemake@params$locus   %||% character(0)
gene_names     <- snakemake@params$gene    %||% character(0)
protein_names  <- snakemake@params$protein %||% character(0)
variant_labels <- snakemake@params$variant_labels


#----- Read each variant CSV, tagging it with the label provided by Snakemake -----#
read_one <- function(path, label) {
  read_csv(path, show_col_types = FALSE) |> mutate(Calibration = label)
}

merged <- map2(snakemake@input$annotated, variant_labels, read_one) |>
  bind_rows() |>
  mutate(Calibration = factor(Calibration, levels = c("uncalibrated", "calibrated")))


#----- Validated-hits filter: each entry is "<target>" (any miRNA) or "<target>,<miRNA>" (only that miRNA) -----#
parse_pairs <- function(entries, target_col) {
  if (length(entries) == 0) {
    return(tibble(!!target_col := character(0), miRNA_filter = character(0)))
  }
  parts <- strsplit(as.character(entries), ",", fixed = TRUE)
  if (any(lengths(parts) > 2L)) {
    bad <- entries[lengths(parts) > 2L][1L]
    stop(sprintf("build_plots %s: expected '<target>' or '<target>,<miRNA>', got: %s", target_col, bad))
  }
  target <- trimws(vapply(parts, `[[`, character(1), 1L))
  mirna  <- trimws(vapply(parts, function(p) if (length(p) > 1L) p[[2L]] else NA_character_, character(1)))
  if (any(target == ""))                stop(sprintf("build_plots %s: empty target before comma in entry: %s", target_col, entries[target == ""][1L]))
  if (any(!is.na(mirna) & mirna == "")) stop(sprintf("build_plots %s: empty miRNA after comma in entry: %s",  target_col, entries[!is.na(mirna) & mirna == ""][1L]))
  tibble(!!target_col := target, miRNA_filter = mirna)
}

locus_pairs   <- parse_pairs(locus_tags,    "locus_tag")
gene_pairs    <- parse_pairs(gene_names,    "gene_name")
protein_pairs <- parse_pairs(protein_names, "protein_name")

if (nrow(protein_pairs) > 0L && !"protein_name" %in% names(merged)) {
  stop("build_plots protein: input CSV has no `protein_name` column.")
}

join_filter <- function(df, pairs, by_col) {
  if (nrow(pairs) == 0L) return(df[0, , drop = FALSE])
  df |>
    inner_join(pairs, by = by_col, relationship = "many-to-many") |>
    filter(is.na(miRNA_filter) | miRNA == miRNA_filter) |>
    select(-miRNA_filter)
}

validated <- bind_rows(
  join_filter(merged, locus_pairs,   "locus_tag"),
  join_filter(merged, gene_pairs,    "gene_name"),
  join_filter(merged, protein_pairs, "protein_name")
)

# One row per (miRNA, locus_tag, Calibration) for the p-value distribution vlines
# and for the legend/scale keys. Position-based plots use the full `validated` set.
validated_best <- validated |>
  group_by(Calibration) |>
  arrange(P_value) |>
  distinct(miRNA, locus_tag, .keep_all = TRUE) |>
  ungroup()

#----- Auto color + label maps for the validated hits, keyed by locus_tag -----#
hit_levels <- unique(validated_best$locus_tag)
hit_colors <- setNames(scales::hue_pal()(length(hit_levels)), hit_levels)
hit_label_table <- validated_best |>
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
  dens <- merged |>
    group_by(Calibration) |>
    summarise(ymax = max(density(P_value, na.rm = TRUE)$y), .groups = "drop")
  ycut     <- ceiling(min(dens$ymax) * 3.5)
  ynatural <- ceiling(max(dens$ymax) * 1.05)

  ggplot(merged, aes(x = P_value)) +
    geom_density(color = "black", fill = "#d1d1d1") +
    geom_vline(
      data = validated_best,
      aes(xintercept = P_value, color = locus_tag),
      linetype = "dashed"
    ) +
    scale_color_manual(name = "Validated hits", values = hit_colors, labels = hit_labels) +
    coord_cartesian(ylim = c(0, ycut)) +
    facet_wrap(~Calibration, labeller = facet_labels) +
    labs(
      x = "p-value", y = "Density",
      title = sprintf("Distribution of p-values (uncut Y-axis: %d)", ynatural)
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
  pick <- if (any(merged$Calibration == "calibrated")) "calibrated" else "uncalibrated"
  base_data <- merged |>
    filter(Calibration == pick) |>
    mutate(facet_label = "With/WithOut Calibration")

  ggplot(base_data, aes(x = Position, y = Energy)) +
    geom_point(
      data = base_data |> filter(P_value <= pvalue_cutoff),
      color = "black", size = 0.7
    ) +
    geom_point(
      data = validated |> filter(Calibration == pick),
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
    facet_wrap(~facet_label) +
    theme_bw(base_size = basesize) +
    theme(strip.placement = "outside")
}


#----- Plot registry: name -> builder -----#
plot_builders <- list(
  pvalue_distribution = plot_pvalue_distribution,
  position_pvalue     = plot_position_pvalue,
  position_energy     = plot_position_energy
)


#----- Resolve the requested selection -----#
# "all" -> every plot in registered order; otherwise comma-split + trimmed.
selected <- if (identical(plot_type, "all")) {
  names(plot_builders)
} else {
  trimws(strsplit(plot_type, ",", fixed = TRUE)[[1]])
}

unknown <- setdiff(selected, names(plot_builders))
if (length(unknown) > 0) {
  stop(sprintf(
    "Unknown build_plots.type entries: %s. Valid: %s, all, or a comma-separated subset.",
    paste(unknown, collapse = ", "),
    paste(names(plot_builders), collapse = ", ")
  ))
}

plots <- lapply(selected, function(name) plot_builders[[name]]())


#----- Compose & save -----#
# Page width: 9.5 in (single), 16 in (double), 18 in (triple = legacy "all").
# Layout widths: equal for the 2-plot case; (2, 2, 1) for the 3-plot case to
# preserve the legacy arranger.R aesthetic where position_energy is narrower.
n <- length(plots)
if (n == 1) {
  combined  <- plots[[1]]
  out_width <- 9.5
} else {
  layout_widths <- if (n == 3) c(2, 2, 1) else rep(1, n)
  out_width     <- if (n == 2) 16 else 18
  combined      <- (Reduce(`|`, plots) + plot_layout(widths = layout_widths, guides = "collect")) & theme(legend.position = "bottom")
}

dir.create(dirname(snakemake@output$pdf), showWarnings = FALSE, recursive = TRUE)
ggsave(
  snakemake@output$pdf, combined,
  width = out_width, height = 5, units = "in"
)