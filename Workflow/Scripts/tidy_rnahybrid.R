library(tidyverse)


parse_args <- function(args) {
  parsed <- list(input = NULL, output = NULL)

  i <- 1
  while (i <= length(args)) {
    arg <- args[[i]]

    if (arg %in% c("-i", "--input")) {
      i <- i + 1
      parsed$input <- args[[i]]
    } else if (arg %in% c("-o", "--output")) {
      i <- i + 1
      parsed$output <- args[[i]]
    } else {
      stop(sprintf("Unknown argument: %s", arg), call. = FALSE)
    }
    i <- i + 1
  }

  if (is.null(parsed$input) || is.null(parsed$output)) {
    stop("Both --input and --output are required.", call. = FALSE)
  }

  parsed
}


main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))

  columns <- c(
    "Gene",
    "Gene_length",
    "miRNA",
    "miRNA_length",
    "Energy",
    "P_value",
    "Position",
    "Target_unmatches",
    "Target_matches",
    "miRNA_matches",
    "miRNA_unmatches"
  )

  df <- tibble(read.csv2(args$input, header = FALSE, sep = ':', col.names = columns))

  df |>
  mutate(Energy = as.double(Energy), P_value = as.double(P_value)) |>
  arrange(Energy) |>
  select(miRNA, Gene, Energy, P_value, miRNA_unmatches, miRNA_matches, Target_matches, Target_unmatches) -> extracted

  write_csv(extracted, args$output)
}


main()