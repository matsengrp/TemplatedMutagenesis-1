#!/usr/bin/env Rscript

# this script takes two csvs, presumably with the false positive rates
# of a set of sequences when compared with different references and
# makes a plot comparing the rates against those references
library(argparse)
library(reshape2)
library(ggplot2)
parser = ArgumentParser()
parser$add_argument("--input-1", dest = "input_1")
parser$add_argument("--input-2", dest = "input_2")
parser$add_argument("--output", dest = "output")
args = parser$parse_args()

fpr1 = read.csv(args$input_1)[,c("hit_fraction", "input_file", "reference", "k")]
fpr2 = read.csv(args$input_2)[,c("hit_fraction", "input_file", "reference", "k")]
fpr = rbind(fpr1, fpr2)
pdf(args$output, width=6, height=4)
ggplot(fpr) +
    geom_point(aes(x = k, y = hit_fraction, color = reference)) +
    xlab("Minimum Template Size") +
    ylab("Fraction of Mutations with Templates") +
    labs(color = "Reference\nSet")
dev.off()
