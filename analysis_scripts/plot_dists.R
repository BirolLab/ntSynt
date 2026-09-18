library(ggplot2)
library(dplyr)
library(tidyr)
library(argparse)
library(ggiraph)
library(htmlwidgets)

# ---------------------------------------------------------------------------
# Load data
# Expected to be run from the mash/ directory (Snakemake cd's there).
# ---------------------------------------------------------------------------

parser <- ArgumentParser(description = "Plot the pairwise Mash distances")
parser$add_argument("--full", help = "Full assembly Mash distances", required = TRUE)
parser$add_argument("--syntenic", help="Syntenic regions Mash distances", required = TRUE)
parser$add_argument("--non-syntenic", help = "Non-syntenic regions Mash distances", required = TRUE)
parser$add_argument("--name-conversions", help = "Name conversion TSV", required = TRUE)
parser$add_argument("-o", help = "Output file prefix", required = TRUE)

args <- parser$parse_args()

col_names <- c("query", "reference", "mash_dist", "p_value", "matching_hashes")

# Read a Mash TSV, returning an empty (correctly-typed) data frame if the file is empty
read_mash_tsv <- function(path) {
  if (file.info(path)$size == 0) {
    return(setNames(
      data.frame(character(), character(), numeric(), numeric(), character(),
                 stringsAsFactors = FALSE),
      col_names
    ))
  }
  read.table(path, col.names = col_names)
}

full_genome  <- read_mash_tsv(args$full)
syntenic     <- read_mash_tsv(args$syntenic)
non_syntenic <- read_mash_tsv(args$non_syntenic)

# Name conversion TSV: old basename -> new (human-readable) genome name.
# Adjust the path below if the file lives somewhere other than the mash/ dir.
name_conversion_file <- args$name_conversions
 
name_conversion <- read.table(name_conversion_file, sep = "\t",
                               col.names = c("old", "new"),
                               stringsAsFactors = FALSE)
 
# Underscores -> spaces for display purposes
name_map <- setNames(gsub("_", " ", name_conversion$new), name_conversion$old)
 
# Translate a vector of cleaned basenames to their display names.
# Anything not found in the lookup is left as-is (with a warning) so the
# script doesn't silently drop data if the conversion table is incomplete.
translate_name <- function(x) {
  if (length(x) == 0) return(character(0))
  translated <- unname(name_map[x])
  missing <- x[is.na(translated)]
  if (length(missing) > 0) {
    warning(sprintf(
      "No name conversion found for: %s (leaving basename as-is)",
      paste(unique(missing), collapse = ", ")
    ))
  }
  ifelse(is.na(translated), x, translated)
}

# Helper: strip path and suffixes to get a clean genome accession ID
clean_name <- function(x) {
  x <- basename(x)
  x <- sub("\\.(syntenic|non_syntenic)\\.\\S*$", "", x)
  x
}

for (df_name in c("full_genome", "syntenic", "non_syntenic")) {
  df <- get(df_name)
  df$query     <- translate_name(clean_name(df$query))
  df$reference <- translate_name(clean_name(df$reference))
  assign(df_name, df)
}

# Remove self-comparisons
full_genome  <- full_genome[full_genome$query != full_genome$reference, ]
syntenic     <- syntenic[syntenic$query != syntenic$reference, ]
non_syntenic <- non_syntenic[non_syntenic$query != non_syntenic$reference, ]

# Create a pair ID (sorted so order doesn't matter), then deduplicate
make_pair_id <- function(df) {
  df$pair_id <- apply(df[, c("query", "reference")], 1, function(x) {
    paste(sort(x), collapse = " vs ")
  })
  df[!duplicated(df$pair_id), ]
}

full_genome  <- make_pair_id(full_genome)
syntenic     <- make_pair_id(syntenic)
non_syntenic <- make_pair_id(non_syntenic)

# ---------------------------------------------------------------------------
# Combine into a single long data frame
# ---------------------------------------------------------------------------

full_genome$region  <- rep("Full genome", nrow(full_genome))
syntenic$region     <- rep("Syntenic", nrow(syntenic))
non_syntenic$region <- rep("Non-syntenic", nrow(non_syntenic))

combined <- bind_rows(full_genome, syntenic, non_syntenic) %>%
  select(pair_id, region, mash_dist)

combined$region <- factor(combined$region,
                          levels = c("Full genome", "Syntenic", "Non-syntenic"))

# ---------------------------------------------------------------------------
# Paired Wilcoxon signed-rank tests (syntenic vs non-syntenic)
# ---------------------------------------------------------------------------

paired <- bind_rows(syntenic, non_syntenic) %>%
  select(pair_id, region, mash_dist) %>%
  pivot_wider(names_from = region, values_from = mash_dist)

# Ensure both columns exist even if one input was empty
for (col in c("Syntenic", "Non-syntenic")) {
  if (!col %in% names(paired)) paired[[col]] <- NA_real_
}

# Only pairs with both measurements present can be used
complete_pairs <- paired %>%
  filter(!is.na(Syntenic), !is.na(`Non-syntenic`))

if (nrow(complete_pairs) >= 3) {
  test_result <- tryCatch(
    wilcox.test(
      complete_pairs$Syntenic,
      complete_pairs$`Non-syntenic`,
      paired      = TRUE,
      exact       = FALSE,
      alternative = "two.sided"
    ),
    error = function(e) {
      warning(sprintf("Wilcoxon test failed (%s); skipping.", conditionMessage(e)))
      NULL
    }
  )

  if (!is.null(test_result)) {
    cat("Paired Wilcoxon signed-rank test (Syntenic vs Non-syntenic)\n")
    cat(sprintf("  V = %.4g,  p-value = %.4g\n\n",
                test_result$statistic, test_result$p.value))

    p_label <- ifelse(
      test_result$p.value < 0.001,
      sprintf("p = %.2e", test_result$p.value),
      sprintf("p = %.4f", test_result$p.value)
    )
  } else {
    p_label <- NULL
  }
} else {
  warning(sprintf(
    "Only %d complete pairs available (need >= 3); skipping Wilcoxon test.",
    nrow(complete_pairs)
  ))
  p_label <- NULL
}

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------

fill_vals   <- c("Full genome"  = "#74C476",
                 "Syntenic"     = "#4393C3",
                 "Non-syntenic" = "#D6604D")
colour_vals <- c("Full genome"  = "#238B45",
                 "Syntenic"     = "#2166AC",
                 "Non-syntenic" = "#B2182B")

# Tooltip shown on hover: which two genomes this point represents
combined$tooltip_label <- paste0(combined$pair_id,
                                  "\nMash distance: ", round(combined$mash_dist * 100, 2), "%")

p <- ggplot(combined, aes(x = region, y = mash_dist*100, fill = region)) +

  # Boxplot (outliers hidden — raw data shown below)
  geom_boxplot(outlier.shape = NA, alpha = 0.6, width = 0.45) +

  # Individual observations
  geom_jitter_interactive(
    aes(colour = region, tooltip = tooltip_label, data_id = pair_id),
    position = position_jitter(width = 0.08, height = 0, seed = 42),
    size = 2.5, alpha = 0.3
  ) +

  scale_fill_manual(values   = fill_vals) +
  scale_colour_manual(values = colour_vals) +

  labs(x = NULL, y = "Mash distance (%)") +

  theme_bw(base_size = 18) +
  theme(legend.position = "none")

# Annotate Wilcoxon p-value between Syntenic and Non-syntenic bars (positions 2 & 3)
if (!is.null(p_label)) {
  y_ann <- max(combined$mash_dist, na.rm = TRUE) * 1.05
  p <- p + annotate("text", x = 2.5, y = y_ann, label = p_label, size = 4)
}

ggsave(paste0(args$o, ".pdf"), p, width = 6, height = 5)
ggsave(paste0(args$o, ".png"), p, width = 6, height = 5, dpi = 300)

girafe_obj <- girafe(
  ggobj = p,
  width_svg = 6,
  height_svg = 5,
  options = list(
    opts_hover(css = "stroke:black;stroke-width:2px;opacity:1;"),
    opts_hover_inv(css = "opacity:0.2;"),
    opts_tooltip(css = paste(
      "background-color:white;",
      "padding:5px 8px;",
      "border-radius:4px;",
      "border:1px solid #888;",
      "font-family:sans-serif;",
      "font-size:12px;",
      "white-space:pre-line;"
    )),
    opts_sizing(rescale = FALSE)
  )
)
 
saveWidget(girafe_obj, paste0(args$o, ".html"), selfcontained = TRUE)

cat("Plots saved to mash_divergence_boxplot.pdf / .png / .html\n")
