library(ggplot2)
reference_name_map = list()
reference_name_map[["mus_musculus_129S1_v_genes"]] = "Mouse V Genes"
reference_name_map[["gpt_132"]] = "gpt Genes"
tissue_types = read.csv("analysis/tissue_annotations.csv")
rownames(tissue_types) = tissue_types$name
paper_theme = theme(legend.key.width = unit(.2, "cm"),
                    legend.key.height = unit(.75, "cm"),
                    axis.text = element_text(size = 14),
                    axis.title = element_text(size = 14),
                    legend.title = element_text(size = 14),
                    legend.text = element_text(size = 12))


## adapted from phyloseq
plot_tree = function(phytree,
                     label.tips=NULL, text.size=NULL,
                     sizebase=5, base.spacing = 0.02,
                     ladderize=FALSE, plot.margin=0.2, title=NULL,
                     treetheme=NULL, justify="jagged"){
    method = "treeonly"
    treeSegs = tree_layout(phytree, ladderize=ladderize)
    edgeMap = aes(x=xleft, xend=xright, y=y, yend=y)
    vertMap = aes(x=x, xend=x, y=vmin, yend=vmax)
    ## Initialize phylogenetic tree.
    ## Naked, lines-only, unannotated tree as first layers. Edge (horiz) first, then vertical.
    p = ggplot(data=treeSegs$edgeDT) + geom_segment(edgeMap) +
        geom_segment(vertMap, data=treeSegs$vertDT)
    ## If no text.size given, calculate it from number of tips ("species", aka taxa)
    ## This is very fast. No need to worry about whether text is printed or not.
    if(is.null(text.size)){
        text.size = manytextsize(ape::Ntip(phytree))
    }
    ## Add the species labels to the right.
    if(!is.null(label.tips) & method!="sampledodge") {
        ## If method is sampledodge, then labels are added to the right of points, later.
        ## Add labels layer to plotting object.
        labelDT = treeSegs$edgeDT[!is.na(OTU), ]
        if(justify == "jagged") {
            ## Tip label aesthetic mapping.
            ## Aesthetics can be NULL, and that aesthetic gets ignored.
            labelMap = aes_string(x="xright", y="y", label=label.tips, color=color)
        } else {
            ## The left-justified version of tip-labels.
            labelMap = aes_string(x="max(xright, na.rm=TRUE)", y="y", label=label.tips, color=color)
        }
        p = p + geom_text(labelMap, data=labelDT, size=I(text.size), hjust=-0.1, na.rm=TRUE)
    }
    ## Theme specification
    if(is.null(treetheme)) {
        ## If NULL, then use the default tree theme.
        treetheme = theme(axis.ticks = element_blank(),
                          axis.title.x=element_blank(), axis.text.x=element_blank(),
                          axis.title.y=element_blank(), axis.text.y=element_blank(),
                          panel.background = element_blank(),
                          panel.grid.minor = element_blank(),
                          panel.grid.major = element_blank())
    }
    if(inherits(treetheme, "theme")) {
        ## If a theme, add theme layer to plot.
        ## For all other cases, skip this, which will cause default theme to be used
        p = p + treetheme
    }
    ## Optionally add a title to the plot
    if(!is.null(title)) {
        p = p + ggtitle(title)
    }
    return(p)
}

tree_layout = function(phy, ladderize=FALSE){
    if(is.null(phy$edge.length)){
        ## If no edge lengths, set them all to value of 1 (dendrogram).
        phy$edge.length = rep(1L, times=nrow(phy$edge))
    }
    ## Perform ladderizing, if requested
    if(ladderize != FALSE) {
        if(ladderize == "left") {
            phy = ape::ladderize(phy, FALSE)
        } else if(ladderize==TRUE | ladderize=="right") {
            phy = ape::ladderize(phy, TRUE)
        } else {
            stop("You did not specify a supported option for argument `ladderize`.")
        }
    }
    ## 'z' is the tree in postorder order used in calls to .C
    ## Descending order of left-hand side of edge (the ancestor to the node)
    z = ape::reorder.phylo(phy, order="postorder")
    ## Initialize some characteristics of the tree.
    Ntip = length(phy$tip.label)
    ROOT = Ntip + 1
    nodelabels = phy$node.label
    ## Horizontal positions
    xx = ape::node.depth.edgelength(phy)
    ## vertical positions
    yy = ape::node.height(phy = phy, clado.style = FALSE)
    ## Initialize an edge data.table
    ## Don't set key, order matters
    edgeDT = data.table::data.table(phy$edge, edge.length=phy$edge.length, OTU=NA_character_)
    ## Add tip.labels if present
    if(!is.null(phy$tip.label)) {
        ## Initialize OTU, set node (V2) as key, assign taxa_names as OTU label
        edgeDT[, OTU:=NA_character_]
        data.table::setkey(edgeDT, V2)
        edgeDT[V2 <= Ntip, OTU:=phy$tip.label]
    }
    ## Add the mapping for each edge defined in `xx` and `yy`
    edgeDT[, xleft:=xx[V1]]
    edgeDT[, xright:=xx[V2]]
    edgeDT[, y:=yy[V2]]
    ## Next define vertical segments
    vertDT = edgeDT[, list(x=xleft[1], vmin=min(y), vmax=max(y)), by=V1, mult="last"]
    if(!is.null(phy$node.label)) {
        ## Add non-root node labels to edgeDT
        edgeDT[V2 > ROOT, x:=xright]
        edgeDT[V2 > ROOT, label:=phy$node.label[-1]]
        ## Add root label (first node label) to vertDT
        data.table::setkey(vertDT, V1)
        vertDT[J(ROOT), y:=mean(c(vmin, vmax))]
        vertDT[J(ROOT), label:=phy$node.label[1]]
    }
    return(list(edgeDT=edgeDT, vertDT=vertDT))
}

manytextsize = function(n, mins=0.5, maxs=4, B=6, D=100){
    ## empirically selected size-value calculator.
    s = B * exp(-n/D)
    ## enforce a floor.
    s = ifelse(s > mins, s, mins)
    ## enforce a max
    s = ifelse(s < maxs, s, maxs)
    return(s)
}

#' @param tab The content of the table, everything that goes below the
#'     bar.
#' @param colnames The column names, everything that goes above the
#'     bar.
#' @param alignment How the columns should be aligned, so something
#'     like "l | l l l" would correspond to four columns, a bar
#'     between the first and second column, and everything left
#'     aligned.
format_table <- function(tab, colnames, alignment) {
    header = paste("\\begin{tabular}{", alignment, "}", sep = "")
    topline = paste(paste(colnames, collapse = "&"), "\\\\\\hline", sep = "")
    rows = apply(tab, 1, function(x) paste(x, collapse = "&"))
    all_rows = paste(rows, collapse = "\\\\\n")
    footer = "\\end{tabular}"
    out = paste(header, topline, all_rows, footer, sep = "\n")
    return(out)
}
