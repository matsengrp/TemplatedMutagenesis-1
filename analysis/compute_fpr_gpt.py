#!/usr/bin/python

from Bio import SeqIO
from gcgcgc.motif_finder import make_kmer_dictionary
from gcgcgc.motif_finder import indexed_motif_finder
from gcgcgc.motif_finder import hit_fraction
from gcgcgc.process_partis import process_partis_poly
from gcgcgc.process_partis import process_partis
import os
import os.path
import sys
import getopt
import pandas as pd

# parse input arguments
data_directory = ""
output_csv = ""
output_csv_poly = ""
options, remainder = getopt.getopt(sys.argv[1:], "i:j:k:")
for (opt, arg) in options:
    if opt == "-i":
        data_directory = arg
    if opt == "-j":
        output_csv = arg
    if opt == "-k":
        output_csv_poly = arg
# partis annotation files
partis_files = os.listdir(data_directory)
# set of homologous gpt genes to be used as potential conversion templates
reference_gpt = "./data/reference_sets/gpt_132.fasta"
refs_gpt = [r for r in SeqIO.parse(reference_gpt, "fasta")]
# set of mouse V genes to be used as potential conversion templates
reference_v = "./data/reference_sets/mus_musculus_129S1_v_genes.fasta"
refs_v = [r for r in SeqIO.parse(reference_v, "fasta")]


def make_row(hits, k, input_file, reference):
    row = {}
    row["hit_fraction"] = hits
    row["k"] = k
    row["input_file"] = input_file
    row["reference"] = reference
    return(row)


# make k-mer dictionaries from the gpt and v gene references for a range of k
k_list = range(8, 15)
kmer_dicts_gpt = [make_kmer_dictionary(refs_gpt, k) for k in k_list]
kmer_dicts_v = [make_kmer_dictionary(refs_v, k) for k in k_list]
# to store the output of indexed_motif_finder
row_list = []
row_list_poly = []
for p in partis_files:
    mutations = process_partis(os.path.join(data_directory, p))
    mutations_poly = process_partis_poly(os.path.join(data_directory, p),
                                         max_spacing=8, max_mutation_rate=.3)
    imf_out_gpt = [indexed_motif_finder(mutations, kmer_dict, k) for
                   (kmer_dict, k) in zip(kmer_dicts_gpt, k_list)]
    imf_out_v = [indexed_motif_finder(mutations, kmer_dict, k) for
                 (kmer_dict, k) in zip(kmer_dicts_v, k_list)]
    pmf_out_gpt = [indexed_motif_finder(mutations_poly, kmer_dict, k) for
                   (kmer_dict, k) in zip(kmer_dicts_gpt, k_list)]
    pmf_out_v = [indexed_motif_finder(mutations_poly, kmer_dict, k) for
                 (kmer_dict, k) in zip(kmer_dicts_v, k_list)]
    # store the fraction of hits to gpt, hits to v from single or double mutations
    for (gpt, v, gpt_poly, v_poly, k) in zip(imf_out_gpt,
                                             imf_out_v,
                                             pmf_out_gpt,
                                             pmf_out_v,
                                             k_list):
        row_list.append(make_row(hit_fraction(gpt), k, p, "gpt"))
        row_list.append(make_row(hit_fraction(v), k, p, "v"))
        row_list_poly.append(make_row(hit_fraction(gpt_poly), k, p, "gpt"))
        row_list_poly.append(make_row(hit_fraction(v_poly), k, p, "v"))

df = pd.DataFrame(row_list)
df_poly = pd.DataFrame(row_list_poly)
df.to_csv(output_csv)
df_poly.to_csv(output_csv_poly)
