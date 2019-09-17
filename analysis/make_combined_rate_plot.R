#!/usr/bin/env Rscript

# this script takes two csvs, presumably with the false positive rates
# of a set of sequences when compared with different references and
# makes a plot comparing the rates against those references
set.seed(1)
library(argparse)
library(reshape2)
library(ggplot2)
library(magrittr)
library(dplyr)
library(cowplot)
source("analysis/plotting_extras.R")
parser = ArgumentParser()
parser$add_argument("--input-1", dest = "input_1")
parser$add_argument("--input-2", dest = "input_2")
parser$add_argument("--output-file", dest = "output_file")
parser$add_argument("--rc", dest = "rc", default = "False")
args = parser$parse_args()

fpr_1 = read.csv(args$input_1)[,c("hit_fraction", "input_file", "reference", "k",
                                                "reverse_complement", "type", "dale_method")]
fpr_2 = read.csv(args$input_2)[,c("hit_fraction", "input_file", "reference", "k",
                                               "reverse_complement", "type", "dale_method")]
fpr = rbind(fpr_1, fpr_2)
# tissue_types was loaded in with the source("analysis/plotting_extras.R") command
fpr$tissue_type = tissue_types[fpr$input_file,]$tissue
fpr = subset(fpr, dale_method == "False")
# make the confidence intervals
cis = fpr %>%
    group_by_at(c("k", "reference", "reverse_complement", "type")) %>%
    summarise(mean = mean(hit_fraction), se = sd(hit_fraction) / sqrt(length(hit_fraction) - 1))

pdf(args$output_file, width=7.5, height=3.5)
ggplot(subset(fpr, type == "pmf" & reverse_complement == args$rc)) +
    geom_point(aes(x = k, y = hit_fraction, shape = tissue_type, color = reference),
               position = position_dodge(width=.4), alpha = .5, size = 1) +
    xlab("Minimum Donor Tract Size") +
    ylab("Fraction of Mutations\nIdentified as Templated by PyPMF\n") +
    labs(color = "Donor Set", shape = "Tissue") +
    scale_color_discrete(breaks = names(reference_name_map),
                         labels = reference_name_map) +
    scale_x_continuous(breaks = seq(8, 14)) +
    scale_shape_manual("", values = c(24,25)) +
    paper_theme +
    geom_errorbar(aes(x = k, color = reference, ymin = mean - 2 * se, ymax = mean + 2 * se),
                  data = subset(cis, reverse_complement == args$rc & type == "pmf"),
                  position = position_dodge(width=.4), width=.3) +
    geom_point(aes(x = k, y = mean, color = reference, ),
               data = subset(cis, reverse_complement == args$rc & type == "pmf"),
               position = position_dodge(width=.4), size=2.4)
dev.off()
