reference_name_map = list()
reference_name_map[["mus_musculus_129S1_v_genes"]] = "Mouse V\nGenes"
reference_name_map[["gpt_132"]] = "gpt\nGenes"
tissue_types = read.csv("analysis/tissue_annotations.csv")
rownames(tissue_types) = tissue_types$name
