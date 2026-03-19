# This script should be run by a python one with the sole purpose to manage the tibble from RNAhybrid compact output file

library(tidyverse)

search <- c(
  miRNA = 'hsa-miR-1224-5p',
  file_name = './miR_1224_5p.txt'
)

df <- tibble(read.csv2('./output_compact.tsv', header = FALSE, sep = ':', col.names = c('Gene', 'Gene_length', 'miRNA', 'miRNA_length', 'Energy', 'P_value', 'Position', 'miRNA_unmatches', 'miRNA_matches', 'target_matches', 'target_unmatches')))

df |>
  mutate(Energy = as.double(Energy), P_value = as.double(P_value)) |>
  filter(miRNA == search['miRNA'], P_value <= 0.05) |>
  arrange(P_value) |>
  # select(Gene, miRNA, Energy, P_value) |>
  # relocate(miRNA) -> miR_1226_5p
  distinct(Gene, miRNA) |>
  select(Gene) -> miRNA_genes_extracted

# print(miR_1226_5p, n = 1000)

write_csv(miRNA_genes_extracted, search['file_name'], col_names = FALSE)