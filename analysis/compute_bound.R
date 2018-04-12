#!/usr/bin/env Rscript

# this script takes two csvs, presumably with the false positive rates
# of a set of sequences when compared with different references and
# makes a plot comparing the rates against those references
set.seed(1)
library(argparse)
library(reshape2)
library(dplyr)
source("analysis/plotting_extras.R")
parser = ArgumentParser()
parser$add_argument("--ebola-rate", dest = "ebola")
parser$add_argument("--gpt-fpr", dest = "gpt")
parser$add_argument("--output", dest = "output")
args = parser$parse_args()

ebola = read.csv(args$ebola)[,c("hit_fraction", "input_file", "reference", "k")]
gpt = read.csv(args$gpt)[,c("hit_fraction", "input_file", "reference", "k")]
gpt = gpt %>% group_by(k, reference) %>% summarise(hit_fraction = mean(hit_fraction))
out = merge(data.frame(gpt[,c("k", "hit_fraction")]), ebola[,c("k", "hit_fraction")], by = "k", suffixes = c("_gpt", "_ebola"))
out$bound = 1 - (1 - out$hit_fraction_ebola) / (1 - out$hit_fraction_gpt)
write.csv(out, file = args$output)


