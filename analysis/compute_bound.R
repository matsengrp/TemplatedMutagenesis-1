#!/usr/bin/env Rscript

set.seed(1)
library(argparse)
library(reshape2)
library(magrittr)
library(dplyr)
source("analysis/plotting_extras.R")
parser = ArgumentParser()
parser$add_argument("--ebola-rate", dest = "ebola")
parser$add_argument("--gpt-fpr", dest = "gpt")
parser$add_argument("--output-csv", dest = "output_csv")
parser$add_argument("--output-tex", dest = "output_tex")
args = parser$parse_args()

ebola = read.csv(args$ebola)[,c("hit_fraction", "input_file", "reference", "k", "type", "reverse_complement")]
gpt = read.csv(args$gpt)[,c("hit_fraction", "input_file", "reference", "k", "type", "reverse_complement")]
gpt = gpt %>% group_by(k, reference, type, reverse_complement) %>% summarise(hit_fraction = mean(hit_fraction))
## for now just do mf/no reverse complements
gpt = subset(gpt, type == "mf" & reverse_complement == "False")
ebola = subset(ebola, type == "mf" & reverse_complement == "False")
out = merge(data.frame(gpt[,c("k", "hit_fraction")]), ebola[,c("k", "hit_fraction")], by = "k", suffixes = c("_gpt", "_ebola"))
out$bound = 1 - (1 - out$hit_fraction_ebola) / (1 - out$hit_fraction_gpt)
out$bound99 = 1 - (.99 - out$hit_fraction_ebola) / (.99 - out$hit_fraction_gpt)
out$bound95 = 1 - (.95 - out$hit_fraction_ebola) / (.95 - out$hit_fraction_gpt)
out$bound90 = 1 - (.9 - out$hit_fraction_ebola) / (.9 - out$hit_fraction_gpt)
tex_table = format_table(round(out, digits = 3),
                         colnames = c("$k$", "PyMF rate", "PyMF FPR", "UB (1)", "UB (.99)", "UB (.95)", "UB (.9)"),
                         alignment = "l | l l l l l l")
cat(tex_table, file = args$output_tex)
write.csv(out, file = args$output_csv)


