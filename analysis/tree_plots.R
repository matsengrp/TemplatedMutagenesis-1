library(gridExtra)
library(ape)
library(phytools)
library(argparse)
library(ggplot2)
source("analysis/plotting_extras.R")
parser = ArgumentParser()
parser$add_argument("--tree-output", dest = "tree_output")
args = parser$parse_args()


## read in the trees
gpt_250_tree = read.tree("data/reference_sets/scripts/gpt_250.tree")
#mouse_tree = read.tree("data/reference_sets/scripts/mus_musculus_129S1_v_genes.tree", comment.char = "")
human_tree = read.tree("data/reference_sets/scripts/imgt_human_v.tree")
## from http://blog.phytools.org/2012/01/function-to-get-descendant-node-numbers.html
getDescendants = function(tree, node, curr = NULL) {
    if(is.null(curr)) {
        curr = vector()
    }
    daughters = tree$edge[which(tree$edge[,1] == node), 2]
    curr = c(curr,daughters)
    w = which(daughters >= length(tree$tip))
    if(length(w) > 0) {
        for(i in 1:length(w)) {
            curr = getDescendants(tree, daughters[w[i]], curr)
        }
    }
    return(curr)
}

emax = which.max(gpt_250_tree$edge.length)
## this is the daughter node of the longest edge
node = gpt_250_tree$edge[emax, 2]
## this gets all the descendants of a node
desc = getDescendants(gpt_250_tree, node)
## subset to just the leaves, no internal nodes
desc.tips = desc[desc <= length(gpt_250_tree$tip.label)]
## the names of the sequences in that clade
seqs.to.keep = gpt_250_tree$tip.label[desc.tips]
gpt_132_tree = drop.tip(gpt_250_tree, tip = gpt_250_tree$tip.label[-c(desc.tips)])
gpt_132_rooted = midpoint.root(gpt_132_tree)
#mouse_rooted = midpoint.root(mouse_tree)
human_rooted = midpoint.root(human_tree)
tree_theme = theme(axis.ticks.x = element_line(),
                   axis.text.x = element_text(size = 14),
                   plot.title = element_text(size = 16),
                   panel.border = element_blank(),
                   axis.line = element_line(colour = "black"))
## create the plots of the phylogenetic trees
p1 = plot_tree(gpt_132_rooted, ladderize = TRUE) +
    ggtitle("gpt genes") +
    tree_theme + 
    scale_x_continuous(breaks = seq(0, .9, by = .2), limits = c(0, .9))
#p2 = plot_tree(mouse_rooted, ladderize = TRUE, nodelabf = nodeplotblank) +
#    ggtitle("Mouse strain 129S1 V genes") +
#    tree_theme +
#    scale_x_continuous(breaks = seq(0, .9, by = .2), limits = c(0, .9))
p3 = plot_tree(human_rooted, ladderize = TRUE) +
    ggtitle("Human IMGT V genes") +
    tree_theme +
    scale_x_continuous(breaks = seq(0, .9, by = .2), limits = c(0, .9))
pdf(args$tree_output, width = 7, height = 3.5)
grid.arrange(p1, p3, nrow = 1)
dev.off()
