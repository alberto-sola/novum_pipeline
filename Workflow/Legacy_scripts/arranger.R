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
merged <- rbind(uncalib7267, calib7267, uncalib167a, calib167a, uncalib396e, calib396e)
remove(uncalib7267, calib7267, uncalib167a, calib167a, uncalib396e, calib396e)

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
  geom_vline(data = merged |> filter((locus_tag %in% c('LGG_RS02140', 'LGG_RS02130') & miRNA == 'ath-miR167a-5p') |
                                     (locus_tag == 'LGG_RS05490' & miRNA == 'mdo-miR-7267-3p') |
                                     (locus_tag == 'LGG_RS03370' & miRNA == 'gma-miR396e')) |> group_by(Calibration) |> arrange(P_value) |> distinct(miRNA, locus_tag, .keep_all = TRUE),
             aes(xintercept = P_value, color = miRNA),
             linetype = "dashed"
  ) +
  scale_color_manual(name = "Validated hits", values = c("mdo-miR-7267-3p" = "red", "ath-miR167a-5p" = "blue", "gma-miR396e" = "green"),
                                              labels = c("mdo-miR-7267-3p" = "miR-7267 | ycnE", "ath-miR167a-5p" = "miR167a | spaC", "gma-miR396e" = "miR396e | lexA")) +
  coord_cartesian(ylim = c(0, 26)) +
  facet_wrap(~Calibration, labeller = my_labels) +
  labs(x = 'p-value', y = 'Density', title = 'Distribution of p-values (Y-axis ~60)') +
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
  geom_point(data = merged |> filter((locus_tag %in% c('LGG_RS02140', 'LGG_RS02130') & miRNA == 'ath-miR167a-5p') |
                                     (locus_tag == 'LGG_RS05490' & miRNA == 'mdo-miR-7267-3p') |
                                     (locus_tag == 'LGG_RS03370' & miRNA == 'gma-miR396e')) |> group_by(Calibration),
             aes(shape = miRNA),
             color = "red"
  ) +
  scale_shape_manual(name = "Validated hits", values = c("mdo-miR-7267-3p" = 16, "ath-miR167a-5p" = 17, "gma-miR396e" = 15),
                                              labels = c("mdo-miR-7267-3p" = "miR-7267 | ycnE", "ath-miR167a-5p" = "miR167a | spaC", "gma-miR396e" = "miR396e | lexA")) +
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
  geom_point(data = merged |> filter((locus_tag %in% c('LGG_RS02140', 'LGG_RS02130') & miRNA == 'ath-miR167a-5p') |
                                     (locus_tag == 'LGG_RS05490' & miRNA == 'mdo-miR-7267-3p') |
                                     (locus_tag == 'LGG_RS03370' & miRNA == 'gma-miR396e')),
             aes(shape = miRNA),
             color = "red"
  ) +
  scale_shape_manual(name = "Validated hits", values = c("mdo-miR-7267-3p" = 16, "ath-miR167a-5p" = 17, "gma-miR396e" = 15),
                                              labels = c("mdo-miR-7267-3p" = "miR-7267 | ycnE", "ath-miR167a-5p" = "miR167a | spaC", "gma-miR396e" = "miR396e | lexA")) +
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