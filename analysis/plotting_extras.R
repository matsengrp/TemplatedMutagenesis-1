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
