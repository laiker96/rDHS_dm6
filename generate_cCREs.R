#!/usr/bin/env Rscript

library(optparse)

# Define command-line options
option_list <- list(
  make_option(c("-l", "--left"), type="character", help="Path to dataset_left.bed"),
  make_option(c("-r", "--right"), type="character", help="Path to dataset_right.bed"),
  make_option(c("-c", "--center"), type="character", help="Path to dataset_central.bed"),
  make_option(c("-a", "--atac"), type="character", help="Path to reference.bed"),
  make_option(c("-i", "--histone_ids"), type="character", help="Path to histone_data_ids.txt"),
  make_option(c("-x", "--atac_ids"), type="character", help="Path to contexts.csv"),
  make_option(c("-o", "--output"), type="character", help="Output file for cCREs", default="cCREs.bed"),
  make_option(c("-p", "--probability"), type="double", help="Quantile cutoff", default=0.6)
)

# Parse options
opt <- parse_args(OptionParser(option_list=option_list))

# Load input files
histone_left <- read.delim(opt$left, header=FALSE)
histone_right <- read.delim(opt$right, header=FALSE)
histone_center <- read.delim(opt$center, header=FALSE)
atac <- read.delim(opt$atac, header=FALSE)
histone_ids <- read.csv(opt$histone_ids, header=TRUE)
atac_ids <- read.csv(opt$atac_ids, header=TRUE)

# Combine histone signal
histone <- data.frame(
  V1 = histone_left$V1,
  V2 = histone_left$V2,
  V3 = histone_left$V3,
  V4 = histone_left$V4
)

for (i in 5:ncol(histone_left)) {
  histone[[paste0("V", i)]] <- pmax(histone_left[[i]], histone_right[[i]], histone_center[[i]])
}

# Clean up
histone <- histone[,-c(1:4)]
colnames(histone) <- histone_ids$context
colnames(atac) <- c("chr", "start", "end", "id", atac_ids$context)

# Keep only rows with all histone signals present
valid_rows <- apply(histone, 1, all)
histone <- histone[valid_rows, ]
atac <- atac[valid_rows, ]
histone <- log10(histone)

# Z-score normalization
calc_z_scores <- function(column) {
  (column - mean(column)) / sd(column)
}
histone <- as.data.frame(apply(histone, 2, calc_z_scores))

# Generate cCREs
probs <- opt$probability

cCREs <- data.frame(
  chr = atac$chr,
  start = atac$start,
  end = atac$end,
  id = paste("cCRE_", seq_len(nrow(atac)), sep = "")
)

for (context in histone_ids$context) {
  cCREs[[context]] <- ifelse(
    (histone[[context]] >= tapply(histone[[context]], atac[[context]], quantile, probs = probs)[1]) &
      atac[[context]] == 1,
    1, 0
  )
}


# Filter and export
cCREs <- cCREs[rowSums(cCREs[,-c(1:4)]) > 0, ]
print("Found the following number of cCREs in each context:")
print(colSums(cCREs[, -c(1:4)]))
write.table(cCREs, opt$output, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
