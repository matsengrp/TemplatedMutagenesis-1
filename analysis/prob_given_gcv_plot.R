#!/usr/bin/env Rscript

# this script takes two csvs, presumably with the false positive rates
# of a set of sequences when compared with different references and
# makes a plot comparing the rates against those references
library(argparse)
library(reshape2)
library(ggplot2)
library(dplyr)
library(magrittr)
parser = ArgumentParser()
parser$add_argument("--input-1", dest = "input_1")
parser$add_argument("--input-2", dest = "input_2")
parser$add_argument("--output", dest = "output")
args = parser$parse_args()

prob1 = read.csv(args$input_1)
prob2 = read.csv(args$input_2)
probs = rbind(prob1, prob2)

grouped = probs %>% group_by(source, reference, k) %>% summarise(avgprob = mean(prob, na.rm = TRUE))

pdf(args$output, width=6, height=4)
ggplot(grouped) + geom_point(aes(x = k, y = avgprob, color = reference)) + 
    xlab("Minimum Template Size") +
    ylab("Average Probability of Mutation Given GCV") +
    labs(color = "Reference\nSet")
dev.off()
