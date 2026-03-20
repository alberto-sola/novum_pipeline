parse_args <- function(args) {
  parsed <- list(
    input = NULL,
    output = NULL,
    mirna = NULL,
    pvalue_cutoff = NULL
  )

  i <- 1
  while (i <= length(args)) {
    arg <- args[[i]]

    if (arg %in% c("-i", "--input")) {
      i <- i + 1
      parsed$input <- args[[i]]
    } else if (arg %in% c("-o", "--output")) {
      i <- i + 1
      parsed$output <- args[[i]]
    } else if (arg == "--mirna") {
      i <- i + 1
      parsed$mirna <- args[[i]]
    } else if (arg == "--pvalue-cutoff") {
      i <- i + 1
      parsed$pvalue_cutoff <- as.numeric(args[[i]])
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
    "miRNA_unmatches",
    "miRNA_matches",
    "Target_matches",
    "Target_unmatches"
  )

  df <- read.table(
    args$input,
    sep = ":",
    header = FALSE,
    col.names = columns,
    stringsAsFactors = FALSE,
    quote = ""
  )

  df$Energy <- as.numeric(df$Energy)
  df$P_value <- as.numeric(df$P_value)

  if (!is.null(args$mirna)) {
    df <- df[df$miRNA == args$mirna, , drop = FALSE]
  }

  if (!is.null(args$pvalue_cutoff)) {
    df <- df[df$P_value <= args$pvalue_cutoff, , drop = FALSE]
  }

  output_dir <- dirname(args$output)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  write.table(
    df,
    file = args$output,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )
}


main()
