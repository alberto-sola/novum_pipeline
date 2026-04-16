library(tidyverse)
library(ggpubr)
library(patchwork)

#----- Parses tibbles from files -----#
uncalib <- tibble(read.csv("Data/Results/04_15/escherichia_coli_1226/escherichia_coli_wo_calibration/rnahybrid_annotated.csv", header = TRUE, sep = ',')) |>
  mutate(Calibration = "uncalibrated")
calib <- tibble(read.csv("Data/Results/04_15/escherichia_coli_1226/escherichia_coli_w_calibration/rnahybrid_annotated.csv", header = TRUE, sep = ',')) |>
  mutate(Calibration = "calibrated")
merged <- rbind(uncalib, calib)
remove(uncalib, calib)
merged <- merged |> mutate(Calibration = factor(Calibration, levels = c("uncalibrated", "calibrated")))

#----- Assign pretty labels to the Calibration factor-like variable -----#
my_labels <- as_labeller(c("uncalibrated" = "Without Calibration",
                           "calibrated" = "With Calibration"))
basesize <- 12

#----- p-value distribution -----#
merged |>
  ggplot(aes(x = P_value)) +
  geom_density(color = 'black', fill = '#d1d1d1') +
  geom_vline(data = merged |> filter(gene_name == "yegH") |> group_by(Calibration) |> slice_min(P_value, n = 1),
             aes(xintercept = P_value),
             linetype = "dashed",
             color = "black"
  ) +
  labs(x = 'p-value', y = 'Density', title = 'Distribution of p-values (Y-axis limited)') +
  facet_wrap(~Calibration, labeller = my_labels) +
  coord_cartesian(ylim = c(0, 25)) +
  theme_bw(base_size = basesize) +
  theme(strip.placement = "outside") -> p1

#----- Position vs p-value ‒ correlation -----#
merged |>
  filter(Energy <= -25) |>
  mutate(highlight = gene_name == "yegH") |>
  arrange(highlight) |>
  ggplot(aes(x = Position, y = -log10(P_value))) +
  geom_point(aes(color = highlight), size = 0.3) +
  scale_color_manual(name = "Transcripts", values = c("FALSE" = "black", "TRUE" = "red"), labels = c("FALSE" = "Others", "TRUE" = "Validated")) +
  geom_smooth(method = 'lm', alpha = 0.1, linewidth = 0.5) +
  stat_cor(method = "spearman") +
  labs(x = "Relative Position", y = "-log10(p-value)", title = "Position vs p-values (energies <= -25KCal/Mol)") +
  theme_bw(base_size = basesize) +
  facet_wrap(~Calibration, labeller = my_labels) +
  theme(strip.placement = "outside") -> p2

#----- Position vs Energy ‒ correlation -----#
merged |>
  filter(Calibration == "calibrated", P_value <= 0.05) |>
  mutate(highlight = gene_name == "yegH", facet_label = "With/WithOut Calibration") |>
  arrange(highlight) |>
  ggplot(aes(x = Position, y = Energy)) +
  geom_point(aes(color = highlight), size = 0.3) +
  scale_color_manual(name = "Transcripts", values = c("FALSE" = "black", "TRUE" = "red"), labels = c("FALSE" = "Others", "TRUE" = "Validated")) +
  geom_smooth(method = 'lm', alpha = 0.1, linewidth = 0.5) +
  stat_cor(method = "spearman") +
  labs(x = "Relative Position", y = "Energy (KCal/Mol)", title = "Position vs Energy (p-values <= 0.78)") +
  theme_bw(base_size = basesize) +
  facet_wrap(~facet_label) +
  theme(strip.placement = "outside") -> p3


((p1 | p2 | p3) + plot_layout(widths = c(2, 2, 1), guides = "collect")) & theme(legend.position = "bottom")
#----- Saves the plot as PDF file -----#
ggsave('Data/Results/04_15/escherichia_coli_1226/plots1.pdf', width = 16.5, height = 5, units = "in")

(p1 | p2)
#----- Saves the plot as PDF file -----#
ggsave('Data/Results/04_15/escherichia_coli_1226/plots2.pdf', width = 16.5 , height = 5, units = "in")