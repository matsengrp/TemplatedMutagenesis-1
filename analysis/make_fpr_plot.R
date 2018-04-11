#!/usr/bin/env Rscript

# this script takes two csvs, presumably with the false positive rates
# of a set of sequences when compared with different references and
# makes a plot comparing the rates against those references
set.seed(1)
library(argparse)
library(reshape2)
library(ggplot2)
source("analysis/plotting_extras.R")
parser = ArgumentParser()
parser$add_argument("--input-1", dest = "input_1")
parser$add_argument("--input-2", dest = "input_2")
parser$add_argument("--output", dest = "output")
args = parser$parse_args()

fpr1 = read.csv(args$input_1)[,c("hit_fraction", "input_file", "reference", "k")]
fpr2 = read.csv(args$input_2)[,c("hit_fraction", "input_file", "reference", "k")]
fpr = rbind(fpr1, fpr2)
# tissue_types was loaded in with the source("analysis/plotting_extras.R") command
fpr$tissue_type = tissue_types[fpr$input_file,]$tissue

cis = expand.grid(reference = unique(fpr$reference), k = unique(fpr$k))
cis$mean = rep(NA, nrow(cis))
cis$se = rep(NA, nrow(cis))
for(i in 1:nrow(cis)) {
    data_sub = subset(fpr, k == cis[i, "k"] & reference == cis[i, "reference"])
    cis$mean[i] = mean(data_sub$hit_fraction)
    cis$se[i] = sd(data_sub$hit_fraction) / sqrt(nrow(data_sub) - 1)
}

pdf(args$output, width=7.5, height=3.5)
ggplot(fpr) +
    geom_point(aes(x = k, y = hit_fraction, color = reference, shape = tissue_type),
               position = position_dodge(width=.4), alpha = .5, size = 1) +
    xlab("Minimum Donor Tract Size") +
    ylab("Fraction of Mutations with\nConversion Donors\n") +
    labs(color = "Donor Set", shape = "Tissue") +
    scale_color_discrete(breaks = names(reference_name_map),
                         labels = reference_name_map) +
    scale_x_continuous(breaks = seq(8, 14)) +
    scale_shape_manual(values = c(24,25)) +
    paper_theme +
    geom_errorbar(aes(x = k, ymin = mean - 2 * se, ymax = mean + 2 * se, color = reference),
                  data = cis, position = position_dodge(width=.4), width=.3) +
    geom_point(aes(x = k, y = mean, color = reference),
                  data = cis, position = position_dodge(width=.4), size=2.4)
dev.off()
