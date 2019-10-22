#!/usr/bin/env Rscript

set.seed(1)
library(argparse)
library(reshape2)
library(magrittr)
library(dplyr)
source("analysis/plotting_extras.R")
parser = ArgumentParser()
parser$add_argument("--naive-rate", dest = "naive")
parser$add_argument("--gpt-fpr", dest = "gpt")
parser$add_argument("--output-csv", dest = "output_csv")
parser$add_argument("--output-tex", dest = "output_tex")
parser$add_argument("--rc", dest="rc", default="False")
parser$add_argument("--dale-method", dest="dale", default="False")
args = parser$parse_args()

naive = read.csv(args$naive)[,c("dale_method", "hit_fraction", "input_file", "reference", "k", "type", "reverse_complement")]
gpt = read.csv(args$gpt)[,c("dale_method", "hit_fraction", "input_file", "reference", "k", "type", "reverse_complement")]
gpt = gpt %>%
    group_by(dale_method, k, reference, type, reverse_complement) %>%
    summarise(hit_fraction = mean(hit_fraction))
gpt = subset(gpt, type == "pmf" & reverse_complement == args$rc & dale_method == args$dale)
naive = naive %>% group_by(dale_method, k, reference, type, reverse_complement) %>%
    summarise(hit_fraction = mean(hit_fraction))
naive = subset(naive, type == "pmf" & reverse_complement == args$rc & dale_method == args$dale)
out = merge(naive[,c("k", "hit_fraction")],
            gpt[,c("k", "hit_fraction")], by = "k", suffixes = c("_naive", "_gpt"))
get_bound = function(naive, fpr, tpr) {
    bound = 1 - (tpr - naive) / (tpr - fpr)
    bound = ifelse(bound < 0, 0, bound)
}
out$bound = get_bound(naive = out$hit_fraction_naive,
                      fpr = out$hit_fraction_gpt, tpr = 1)
out$bound99 = get_bound(naive = out$hit_fraction_naive,
                      fpr = out$hit_fraction_gpt, tpr = .99)
out$bound95 = get_bound(naive = out$hit_fraction_naive,
                      fpr = out$hit_fraction_gpt, tpr = .95)
out$bound90 = get_bound(naive = out$hit_fraction_naive,
                      fpr = out$hit_fraction_gpt, tpr = .9)

tex_table = format_table(round(out, digits = 2),
                         colnames = c("$k$", "PyPMF rate", "PyPMF FPR", "UB (1)", "UB (.99)", "UB (.95)", "UB (.9)"),
                         alignment = "l | l l l l l l")
cat(tex_table, file = args$output_tex)
write.csv(out, file = args$output_csv)


