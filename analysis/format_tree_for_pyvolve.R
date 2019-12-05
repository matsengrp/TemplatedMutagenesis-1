library(argparse)
library(ape)
parser = ArgumentParser()
parser$add_argument("--input-tree-prefix", dest = "input_tree_prefix")
args = parser$parse_args()


tr = read.tree(paste0(args$input_tree_prefix, ".tree"))
## We need to make sure that all the node names start with letters
tr$node.label = paste("N", tr$node.label, sep = "")
## We scale up the branch lengths
tr$edge.length = 3 * tr$edge.length
write.tree(tr, paste0(args$input_tree_prefix, "_for_pyvolve.tree"))
