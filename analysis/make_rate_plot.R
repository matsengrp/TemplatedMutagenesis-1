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
parser$add_argument("--input", dest = "input")
parser$add_argument("--output-mf", dest = "output_mf")
parser$add_argument("--output-pmf", dest = "output_pmf")
args = parser$parse_args()

fpr = read.csv(args$input)[,c("hit_fraction", "input_file", "reference", "k",
                              "reverse_complement", "type")]
# tissue_types was loaded in with the source("analysis/plotting_extras.R") command
fpr$tissue_type = tissue_types[fpr$input_file,]$tissue

# make the confidence intervals
cis = fpr %>%
    group_by_at(c("k", "reference", "reverse_complement", "type")) %>%
    summarise(mean = mean(hit_fraction), se = sd(hit_fraction) / sqrt(length(hit_fraction) - 1))

pdf(args$output_mf, width=7.5, height=3.5)
ggplot(subset(fpr, type == "mf" & reverse_complement == "False")) +
    geom_point(aes(x = k, y = hit_fraction, shape = tissue_type),
               position = position_dodge(width=.4), alpha = .5, size = 1) +
    xlab("Minimum Donor Tract Size") +
    ylab("Fraction of Mutations with\nConversion Donors\n") +
    labs(color = "Donor Set", shape = "Tissue") +
    scale_color_discrete(breaks = names(reference_name_map),
                         labels = reference_name_map) +
    scale_x_continuous(breaks = seq(8, 14)) +
    scale_shape_manual(values = c(24,25)) +
    paper_theme +
    geom_errorbar(aes(x = k, ymin = mean - 2 * se, ymax = mean + 2 * se),
                  data = subset(cis, reverse_complement == "False" & type == "mf"),
                  position = position_dodge(width=.4), width=.3) +
    geom_point(aes(x = k, y = mean),
               data = subset(cis, reverse_complement == "False" & type == "mf"),
               position = position_dodge(width=.4), size=2.4)
dev.off()

pdf(args$output_pmf, width=7.5, height=3.5)
ggplot(subset(fpr, type == "pmf" & reverse_complement == "False")) +
    geom_point(aes(x = k, y = hit_fraction, shape = tissue_type),
               position = position_dodge(width=.4), alpha = .5, size = 1) +
    xlab("Minimum Donor Tract Size") +
    ylab("Fraction of Mutations with\nConversion Donors\n") +
    labs(color = "Donor Set", shape = "Tissue") +
    scale_color_discrete(breaks = names(reference_name_map),
                         labels = reference_name_map) +
    scale_x_continuous(breaks = seq(8, 14)) +
    scale_shape_manual(values = c(24,25)) +
    paper_theme +
    geom_errorbar(aes(x = k, ymin = mean - 2 * se, ymax = mean + 2 * se),
                  data = subset(cis, reverse_complement == "False" & type == "pmf"),
                  position = position_dodge(width=.4), width=.3) +
    geom_point(aes(x = k, y = mean),
               data = subset(cis, reverse_complement == "False" & type == "pmf"),
               position = position_dodge(width=.4), size=2.4)
dev.off()
