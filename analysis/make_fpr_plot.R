#!/usr/bin/env Rscript
library(reshape2)
library(ggplot2)
args = commandArgs(trailingOnly = TRUE)

if(length(args) != 2) {
    stop("Two arguments (input csv and output pdf) are required.")
}
input = args[1]
output = args[2]

fpr = read.csv(input)[,c("hit_fraction", "input_file", "reference", "k")]
pdf(output, width=6, height=4)
ggplot(fpr) +
    geom_point(aes(x = k, y = hit_fraction, color = reference)) +
    xlab("Minimum Template Size") +
    ylab("Fraction of Mutations with Templates") +
    labs(color = "Reference\nSet")
dev.off()
