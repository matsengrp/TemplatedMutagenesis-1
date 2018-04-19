#!/usr/bin/env Rscript

# this script takes two csvs, presumably with the false positive rates
# of a set of sequences when compared with different references and
# makes a plot comparing the rates against those references
library(argparse)
library(reshape2)
library(ggplot2)
library(dplyr)
library(magrittr)
library(lme4)
library(broom)
library(cowplot)
source("analysis/plotting_extras.R")
parser = ArgumentParser()
parser$add_argument("--input-1", dest = "input_1")
parser$add_argument("--input-2", dest = "input_2")
parser$add_argument("--output", dest = "output")
args = parser$parse_args()

prob1 = read.csv(args$input_1)
prob2 = read.csv(args$input_2)
probs = rbind(prob1, prob2)

# get the average probability for each sample
grouped = probs %>% group_by(source, reference, k) %>% summarise(avgprob = mean(prob, na.rm = TRUE))
# add tissue type
rownames(tissue_types) = paste(rownames(tissue_types), ".csv", sep = "")
grouped$tissue_type = tissue_types[grouped$source,]$tissue


cis = expand.grid(reference = unique(probs$reference), k = unique(probs$k))
cis$estimates = rep(NA, nrow(cis))
cis$ses = rep(NA, nrow(cis))
for(i in 1:nrow(cis)) {
    ref = cis[i,1]
    kp = cis[i,2]
    t = tryCatch({
        out.lmer = lmer(prob ~ 1 + (1 | source/query_name), data = subset(probs, reference == ref & k == kp))
        cis$estimates[i] = tidy(out.lmer)[1,2]
        cis$ses[i] = tidy(out.lmer)[1,3]
    }, error = function(err) {
        return("try-error")
    })
    # for the mouse v genes with k = 14 there isn't enough data to fit the model above, so we do a smaler one
    if(t == "try-error") {
        out.lmer = lmer(prob ~ 1 + (1 | source), data = subset(probs, reference == ref & k == kp))
        cis$estimates[i] = tidy(out.lmer)[1,2]
        cis$ses[i] = tidy(out.lmer)[1,3]
    }
}

pdf(args$output, width=7.5, height=3.5)
ggplot(grouped) +
    geom_point(aes(x = k, y = avgprob, color = reference, shape = tissue_type),
                size = 1, alpha = .6, position = position_dodge(width=.4)) +
    xlab("Minimum Donor Tract Size") +
    ylab("Average Probability of\nMutation Given GCV\n") +
    ylim(c(0, max(grouped$avgprob))) +
    scale_color_discrete(breaks = names(reference_name_map), labels = reference_name_map) +
    scale_shape_manual(values = c(24, 25)) +
    scale_x_continuous(breaks = seq(8, 14)) +
    labs(color = "Donor Set", shape = "Tissue") +
    paper_theme +
    geom_point(aes(x = k, y = estimates, color = reference),
               data = cis, position = position_dodge(width=.4), size = 2.4) +
    geom_errorbar(aes(x = k, ymin = estimates - 2 * ses, ymax = estimates + 2 * ses, color = reference),
                  data = cis, width = .1, position = position_dodge(width=.4), width=.3)
dev.off()

## tests
## try a paired test (this way ends up dropping a lot of data points, particularly for large k)
probs_merged = merge(prob1, prob2, by = c("source", "query_name", "query_mutation_index", "k"))
probs_merged$gpt_minus_v = probs_merged$prob.x - probs_merged$prob.y
stats = list()
for(kp in unique(probs_merged$k)) {
    t = tryCatch({
        stats[[kp]] = tidy(lmer(gpt_minus_v ~ 1 + (1 | source/query_name), data = subset(probs_merged, k == kp)))[1, "statistic"]
    }, error = function(err) {
        return("try-error")
    })
    if(t == "try-error") {
        stats[[kp]] = tidy(lmer(gpt_minus_v ~ 1 + (1 | source), data = subset(probs_merged, k == kp)))[1, "statistic"]
    }
}
## p-values comparing average probability due to gene conversion from the two reference sets for k in 8 to 14
print("p-values for confidence gpt vs. v comparison")
print(sapply(stats[8:14], function(s) pt(s, df = 11, lower.tail = TRUE) * 2))
