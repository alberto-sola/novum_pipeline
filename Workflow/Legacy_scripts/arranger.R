library(tidyverse)
library(ggpubr)
library(patchwork)

#----- Parses tibbles from files -----#
uncalib <- tibble(read.csv("Data/Results/04_15/001_/fusobacterium_nucleatum_wo_calibration/rnahybrid_annotated.csv", header = TRUE, sep = ',')) |>
  mutate(Calibration = "uncalibrated")
calib <- tibble(read.csv("Data/Results/04_15/001_/fusobacterium_nucleatum_w_calibration/rnahybrid_annotated.csv", header = TRUE, sep = ',')) |>
  mutate(Calibration = "calibrated")
merged <- rbind(uncalib, calib)
remove(uncalib, calib)
merged <- merged |> mutate(Calibration = factor(Calibration, levels = c("uncalibrated", "calibrated")))

#----- Assign pretty labels to the Calibration factor-like variable -----#
my_labels <- as_labeller(c("uncalibrated" = "Without Calibration",
                           "calibrated" = "With Calibration"))

#----- p-value distribution -----#
merged |>
  ggplot(aes(x = P_value)) +
  geom_density(color = 'black', fill = '#d1d1d1', linewidth = 0.4) +
  labs(x = 'p-value', y = 'Density', title = 'Distribution of p-values') +
  theme_bw() +
  facet_wrap(~Calibration, labeller = my_labels) +
  theme(strip.placement = "outside") -> p1

#----- Position vs p-value ‒ correlation -----#
merged |>
  mutate(highlight = protein_name == "16S ribosomal RNA") |>
  arrange(highlight) |>
  ggplot(aes(x = Position, y = -log10(P_value))) +
  geom_point(aes(color = highlight), size = 1) +
  scale_color_manual(name = "Transcripts", values = c("FALSE" = "black", "TRUE" = "red"), labels = c("FALSE" = "Others", "TRUE" = "Validated")) +
  geom_smooth(method = 'lm', alpha = 0.1, linewidth = 0.5) +
  stat_cor(method = "spearman") +
  labs(y = "-log10(p-value)", title = "Position vs p-values") +
  theme_bw() +
  facet_wrap(~Calibration, labeller = my_labels) +
  theme(strip.placement = "outside") -> p2

#----- Position vs Energy ‒ correlation -----#
merged |>
  filter(Calibration == "calibrated") |>
  mutate(highlight = protein_name == "16S ribosomal RNA") |>
  arrange(highlight) |>
  ggplot(aes(x = Position, y = Energy)) +
  geom_point(aes(color = highlight), size = 1) +
  scale_color_manual(name = "Transcripts", values = c("FALSE" = "black", "TRUE" = "red"), labels = c("FALSE" = "Others", "TRUE" = "Validated")) +
  geom_smooth(method = 'lm', alpha = 0.1, linewidth = 0.5) +
  stat_cor(method = "spearman") +
  labs(title = "Position vs Energy") +
  theme_bw() -> p3


# (p1 | p2 | p3) | plot_layout(widths = c(2, 2, 1), guides = "collect")
((p1 | p2 | p3) + plot_layout(widths = c(2, 2, 1), guides = "collect")) & theme(legend.position = "bottom")

#----- Saves the plot as PDF file -----#
ggsave('Data/Results/04_15/001_/plots.pdf', width = 9000, height = 2400, units = "px", dpi = 450)