#!/usr/bin/env Rscript

library(ggplot2)
library(reshape2)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2 || args[1] != "-s") {
  stop("Usage: Rscript plot01.R -s <sample_name>")
}
sample <- args[2]

# Construct the file path for the specific sample
fname <- file.path(paste0(sample, '.call.o.txt'))

# Check if file exists
if (!file.exists(fname)) {
  stop(paste("File not found:", fname))
}

# Check if file is empty
info = file.info(fname)
if (info$size == 0) {
  stop(paste("File is empty:", fname))
}

# Create PDF with sample name
pdf_name <- paste0("converted_ratio_", sample, ".pdf")
pdf(pdf_name, width = 3.7, height = 2.5)

# Read and process the data
stats <- read.delim(fname, header = F)
stats <- as.numeric(unlist(stats[1,]))
stats <- matrix(stats, nrow = 8)
stats <- t(stats[c(2,4,6,8),]/(stats[c(1,3,5,7),] + stats[c(2,4,6,8),]))
rownames(stats) <- colnames(stats) <- c('A','T','C','G')
stats <- as.data.frame(stats)
stats$preceding_base <- rownames(stats)

# Create and print the plot
stats <- melt(stats, variable.name = "following_base", value.name = "converted_ratio")
p1 <- ggplot(data = stats, aes(x=following_base, y = preceding_base)) + 
  geom_tile(aes(fill = converted_ratio)) + 
  geom_text(aes(label=round(converted_ratio, 3))) + 
  theme_classic() + 
  theme(axis.line = element_blank(), plot.title = element_text(size = 5, face = "bold")) + 
  scale_fill_gradient(low ='white',high ='red', limits=c(0,1)) + 
  ggtitle(sample)
print(p1)

dev.off()
