library(tidyverse)
library(ggpubr)
library(patchwork)

#----- Parses tibbles from files -----#
uncalib7267 <- tibble(read.csv("Data/Results/lacticaseibacillus_rhamnosus/miR-7267/lacticaseibacillus_rhamnosus_wo_calibration/rnahybrid_annotated.csv", header = TRUE, sep = ',')) |>
  mutate(Calibration = "uncalibrated")
calib7267 <- tibble(read.csv("Data/Results/lacticaseibacillus_rhamnosus/miR-7267/lacticaseibacillus_rhamnosus_w_calibration/rnahybrid_annotated.csv", header = TRUE, sep = ',')) |>
  mutate(Calibration = "calibrated")
uncalib167a <- tibble(read.csv("Data/Results/lacticaseibacillus_rhamnosus/miR167a/lacticaseibacillus_rhamnosus_wo_calibration/rnahybrid_annotated.csv", header = TRUE, sep = ',')) |>
  mutate(Calibration = "uncalibrated")
calib167a <- tibble(read.csv("Data/Results/lacticaseibacillus_rhamnosus/miR167a/lacticaseibacillus_rhamnosus_w_calibration/rnahybrid_annotated.csv", header = TRUE, sep = ',')) |>
  mutate(Calibration = "calibrated")
uncalib396e <- tibble(read.csv("Data/Results/lacticaseibacillus_rhamnosus/miR396e/lacticaseibacillus_rhamnosus_wo_calibration/rnahybrid_annotated.csv", header = TRUE, sep = ',')) |>
  mutate(Calibration = "uncalibrated")
calib396e <- tibble(read.csv("Data/Results/lacticaseibacillus_rhamnosus/miR396e/lacticaseibacillus_rhamnosus_w_calibration/rnahybrid_annotated.csv", header = TRUE, sep = ',')) |>
  mutate(Calibration = "calibrated")

#----- Merges the tables and purges unused objects -----#
merged <- rbind(uncalib, calib) #uncalib167a, calib167a, uncalib396e, calib396e)
remove(uncalib, calib) #uncalib167a, calib167a, uncalib396e, calib396e)

#----- Data tidying -----#
merged <- merged |> mutate(Calibration = factor(Calibration, levels = c("uncalibrated", "calibrated")))

#----- Assign pretty labels to the Calibration factor-like variable -----#
my_labels <- as_labeller(c("uncalibrated" = "Without Calibration",
                           "calibrated" = "With Calibration"))

#----- Plots' text size constant -----#
basesize <- 12

#----- p-value distribution -----#
merged |>
  ggplot(aes(x = P_value)) +
  geom_density(color = 'black', fill = '#d1d1d1') +
  geom_vline(data = merged |> filter((gene_name %in% c('trpA', 'trpB', 'trpC', 'trpD')) |
                                     (locus_tag %in% c('FE838_RS16085', 'FE838_RS16065', 'FE838_RS16080'))) |> group_by(Calibration) |> arrange(P_value) |> distinct(miRNA, locus_tag, .keep_all = TRUE),
             aes(xintercept = P_value, color = locus_tag),
             linetype = "dashed"
  ) +
  scale_color_manual(name = "Validated hits", values = c("FE838_RS16060" = "black", "FE838_RS16090" = "blue", "FE838_RS16070" = "green", "FE838_RS16075" = "lightblue", "FE838_RS16085" = "orange", "FE838_RS16065" = "red", "FE838_RS16080" = "magenta"),
                     labels = c("FE838_RS16060" = "miR-21 | trpA", "FE838_RS16090" = "miR-21 | trpB", "FE838_RS16070" = "miR-21 | trpC", "FE838_RS16075" = "miR-21 | trpD", "FE838_RS16085" = "miR-21 | trpE", "FE838_RS16065" = "miR-21 | trpF", "FE838_RS16080" = "miR-21 | trpG")) +
  coord_cartesian(ylim = c(0, 50)) +
  facet_wrap(~Calibration, labeller = my_labels) +
  labs(x = 'p-value', y = 'Density', title = 'Distribution of p-values (Y-axis 800+)') +
  theme_bw(base_size = basesize) +
  theme(strip.placement = 'outside') -> p1

#----- Position vs p-value ‒ correlation -----#
merged |>
  ggplot(aes(x = Position, y = -log10(P_value))) +
  geom_point(data = merged |> filter(P_value <= 0.01) |>
    group_by(Calibration),
             # aes(color = "non_validated")
             color = "black",
             size = 0.7
  ) +
  geom_point(data = merged |> filter((gene_name %in% c('trpA', 'trpB', 'trpC', 'trpD')) |
                                     (locus_tag %in% c('FE838_RS16085', 'FE838_RS16065', 'FE838_RS16080'))) |> group_by(Calibration),
             aes(shape = locus_tag),
             color = "red"
  ) +
  scale_shape_manual(name = "Validated hits", values = c("FE838_RS16060" = 16, "FE838_RS16090" = 17, "FE838_RS16070" = 15, "FE838_RS16075" = 19, "FE838_RS16085" = 20, "FE838_RS16065" = 21, "FE838_RS16080" = 22),
                     labels = c("FE838_RS16060" = "miR-21 | trpA", "FE838_RS16090" = "miR-21 | trpB", "FE838_RS16070" = "miR-21 | trpC", "FE838_RS16075" = "miR-21 | trpD", "FE838_RS16085" = "miR-21 | trpE", "FE838_RS16065" = "miR-21 | trpF", "FE838_RS16080" = "miR-21 | trpG")) +
  geom_smooth(method = 'lm', linewidth = 0.5) +
  stat_cor(method = "spearman") +
  labs(x = "Relative Position", y = "-log10(p-value)", title = "Position vs p-values (\U2264 0.01)") +
  theme_bw(base_size = basesize) +
  facet_wrap(~Calibration, labeller = my_labels) +
  theme(strip.placement = "outside") -> p2

#----- Position vs energy -----#
merged |> filter(Calibration == 'calibrated') |> mutate(facet_label = "With/WithOut Calibration") |>
  ggplot(aes(x = Position, y = Energy)) +
  geom_point(data = merged |> filter(P_value <= 0.01) |>
    group_by(Calibration),
             color = 'black',
             size = 0.7
  ) +
  geom_point(data = merged |> filter((gene_name %in% c('trpA', 'trpB', 'trpC', 'trpD')) |
                                     (locus_tag %in% c('FE838_RS16085', 'FE838_RS16065', 'FE838_RS16080'))),
             aes(shape = locus_tag),
             color = "red"
  ) +
  scale_shape_manual(name = "Validated hits", values = c("FE838_RS16060" = 16, "FE838_RS16090" = 17, "FE838_RS16070" = 15, "FE838_RS16075" = 19, "FE838_RS16085" = 20, "FE838_RS16065" = 21, "FE838_RS16080" = 22),
                     labels = c("FE838_RS16060" = "miR-21 | trpA", "FE838_RS16090" = "miR-21 | trpB", "FE838_RS16070" = "miR-21 | trpC", "FE838_RS16075" = "miR-21 | trpD", "FE838_RS16085" = "miR-21 | trpE", "FE838_RS16065" = "miR-21 | trpF", "FE838_RS16080" = "miR-21 | trpG")) +
  geom_smooth(method = 'lm', linewidth = 0.5) +
  stat_cor(method = "spearman") +
  labs(x = "Relative Position", y = "Energy (KCal/Mol)", title = "Position vs Energy (p-value \U2264 0.01)") +
  facet_wrap(~facet_label) +
  theme_bw(base_size = basesize) +
  theme(strip.placement = "outside") -> p3

#----- Saves the plot as PDF file -----#
((p1 | p2 | p3) + plot_layout(widths = c(2, 2, 1), guides = "collect")) & theme(legend.position = "bottom")
ggsave('Data/Results/lacticaseibacillus_rhamnosus/plots1.pdf', width = 18, height = 5, units = "in")

(p1 | p2) + plot_layout(guides = 'collect') & theme(legend.position = "bottom")
ggsave('Data/Results/lacticaseibacillus_rhamnosus/plots2.pdf', width = 18, height = 5, units = "in")