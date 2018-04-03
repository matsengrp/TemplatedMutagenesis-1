#!/usr/bin/env Rscript

# this script takes two csvs, with the per base mutation probabilities
# given gene conversion of a set of sequences when compared with
# different references and makes a plot comparing the rates against
# those references
library(argparse)
library(reshape2)
library(ggplot2)
parser = ArgumentParser()
parser$add_argument("--input-1", dest = "input_1")
parser$add_argument("--input-2", dest = "input_2")
parser$add_argument("--output", dest = "output")
args = parser$parse_args()

df1 = read.csv(args$input_1)
df2 = read.csv(args$input_2)
perbase = rbind(df1, df2)

expected = t(apply(perbase, 1, function(x) {
    gl = x[2]
    mutated = x[3]
    templates = as.numeric(x[8:11])
    names(templates) = names(x)[8:11]
    templates[which(names(templates) == gl)] = NA
    return(templates / sum(templates, na.rm = TRUE))
}))
expected = data.frame(expected)
expected$germline = perbase$gl_base
expected$reference = perbase$reference
expected$k = perbase$k
expected$id = 1:nrow(expected)

observed = t(apply(perbase, 1, function(x) {
    gl = x[2]
    mutated = x[3]
    out = numeric(4)
    names(out) = c("A", "C", "G", "T")
    out[which(names(out) == mutated)] = 1
    return(out)
}))
observed = data.frame(observed)
observed$germline = perbase$gl_base
observed$reference = perbase$reference
observed$k = perbase$k
observed$id = 1:nrow(observed)

observedm = melt(observed, id.vars = c("id", "germline", "reference", "k"))
colnames(observedm)[5:6] = c("target", "observed")
expectedm = melt(expected, id.vars = c("id", "germline", "reference", "k"))
colnames(expectedm)[5:6] = c("target", "probability")
full = merge(observedm, expectedm, by = c("id", "germline", "target", "reference", "k"))
pdf(args$output, width = 8, height = 5)
ggplot(subset(full, reference == "gpt_132" & k == 8), aes(x = probability, y = observed)) +
    geom_jitter(alpha = .01, size = .5, width = .02, height = .02) +
    facet_grid(germline ~ target) +
    stat_smooth(method = "loess", se = FALSE) +
    xlab("Probability of Mutation Under GCV Model") +
    ylab("Mutation Observed") +
    ggtitle("gpt reference, k = 8")

ggplot(subset(full, reference == "mus_musculus_129S1_v_genes" & k == 8), aes(x = probability, y = observed)) +
    geom_jitter(alpha = .01, size = .5, width = .02, height = .02) +
    facet_grid(germline ~ target) +
    stat_smooth(method = "loess", se = FALSE) +
    xlab("Probability of Mutation Under GCV Model") +
    ylab("Mutation Observed") +
    ggtitle("Mouse reference, k = 8")
dev.off()
