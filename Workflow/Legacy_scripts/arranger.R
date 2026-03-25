library(tidyverse)

# var_names <- c(
#   'Gene',
#   'Gene_length',
#   'miRNA',
#   'miRNA_length',
#   'Energy',
#   'P_value',
#   'Position',
#   'miRNA_unmatches',
#   'miRNA_matches',
#   'target_matches',
#   'target_unmatches'
# )

df <- tibble(read.csv("Data/Results/03_25/007_/escherichia_coli/rnahybrid_annotated.csv", header = TRUE, sep = ','))

glimpse(df |>
  mutate(Energy = as.double(Energy), P_value = as.double(P_value)) |>
  arrange(Energy) |>
  select(miRNA, Gene, Energy, P_value)
)
