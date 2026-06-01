library(ggplot2)
library(dplyr)
library(tidyr)

# ---------------------------------------------------------------------------
# Load data
# Expected to be run from the mash/ directory (Snakemake cd's there).
# ---------------------------------------------------------------------------

col_names <- c("query", "reference", "mash_dist", "p_value", "matching_hashes")

full_genome  <- read.table("full_assemblies_dist.tsv", col.names = col_names)
syntenic     <- read.table("synteny_dists.tsv",        col.names = col_names)
non_syntenic <- read.table("non_synteny_dists.tsv",    col.names = col_names)

# Helper: strip path and suffixes to get a clean genome accession ID
clean_name <- function(x) {
  x <- basename(x)
  x <- sub("\\.chr\\.(syntenic|non_syntenic)\\.fa$", "", x)
  x <- sub("_genomic\\.fna$", "", x)
  x
}

for (df_name in c("full_genome", "syntenic", "non_syntenic")) {
  df <- get(df_name)
  df$query     <- clean_name(df$query)
  df$reference <- clean_name(df$reference)
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

full_genome$region  <- "Full genome"
syntenic$region     <- "Syntenic"
non_syntenic$region <- "Non-syntenic"

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

if (nrow(paired) >= 3) {
  test_result <- wilcox.test(
    paired$Syntenic,
    paired$`Non-syntenic`,
    paired      = TRUE,
    exact       = FALSE,
    alternative = "two.sided"
  )
  cat("Paired Wilcoxon signed-rank test (Syntenic vs Non-syntenic)\n")
  cat(sprintf("  V = %.4g,  p-value = %.4g\n\n",
              test_result$statistic, test_result$p.value))

  p_label <- ifelse(
    test_result$p.value < 0.001,
    sprintf("p = %.2e", test_result$p.value),
    sprintf("p = %.4f", test_result$p.value)
  )
} else {
  warning("Too few pairs for Wilcoxon test; skipping.")
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

p <- ggplot(combined, aes(x = region, y = mash_dist, fill = region)) +

  # Boxplot (outliers hidden — raw data shown below)
  geom_boxplot(outlier.shape = NA, alpha = 0.6, width = 0.45) +

  # Individual observations
  geom_jitter(aes(colour = region),
              width = 0.08, height = 0, size = 2.5, alpha = 0.3) +

  scale_fill_manual(values   = fill_vals) +
  scale_colour_manual(values = colour_vals) +

  labs(x = NULL, y = "Mash distance") +

  theme_bw(base_size = 13) +
  theme(legend.position = "none")

# Annotate Wilcoxon p-value between Syntenic and Non-syntenic bars (positions 2 & 3)
if (!is.null(p_label)) {
  y_ann <- max(combined$mash_dist, na.rm = TRUE) * 1.05
  p <- p + annotate("text", x = 2.5, y = y_ann, label = p_label, size = 4)
}

ggsave("mash_divergence_boxplot.pdf", p, width = 6, height = 5)
ggsave("mash_divergence_boxplot.png", p, width = 6, height = 5, dpi = 300)

cat("Plots saved to mash_divergence_boxplot.pdf / .png\n")
