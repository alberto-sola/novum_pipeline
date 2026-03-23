library(tidyverse)



df <- tibble(read.csv2("Data/Results/rnahybrid_output.tsv", header = FALSE, sep = ':', col.names = c(
  'Gene',
  'Gene_length',
  'miRNA',
  'miRNA_length',
  'Energy',
  'P_value',
  'Position',
  'miRNA_unmatches',
  'miRNA_matches',
  'target_matches',
  'target_unmatches')
  )
)

glimpse(df |>
  mutate(Energy = as.double(Energy), P_value = as.double(P_value)) |>
  arrange(Energy) |>
  select(miRNA, Gene, Energy, P_value)
)
