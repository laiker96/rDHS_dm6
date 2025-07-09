#! /usr/bin/env Rscript
#Author:laiker96
suppressMessages(library(optparse))
# Argument parser
option_list = list(
  make_option(c('-f', '--file'), type = 'character'
              ,default = NULL
              ,help = 'Input file name'
              ,metavar = 'character'),
  make_option(c('-o', '--outname'), type = 'character'
              ,default = NULL
              ,help = 'Output file name'
              ,metavar = 'character')
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

df <- read.table(file = opt$file, header = F, sep = "\t")
outname <- opt$outname

process_qpois_values <- function(qpois_values) {
  
  return(10 ** -(qpois_values))

}

df$V4 <- process_qpois_values(df$V4)
write.table(x = df, file = outname, quote = F, sep = "\t", row.names = F, col.names = F)