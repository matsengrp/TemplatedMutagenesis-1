#!/usr/bin/python

from Bio import SeqIO
from gcgcgc.motif_finder import make_kmer_dictionary
from gcgcgc.motif_finder import indexed_motif_finder
from gcgcgc.motif_finder import hit_fraction
from gcgcgc.process_partis import process_partis_poly
from gcgcgc.process_partis import process_partis
import sys
import os
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
partis_files = os.listdir(data_directory)
references = "./data/reference_sets/imgt_ighv_human.fasta"
refs = [r for r in SeqIO.parse(references, "fasta")]


def make_row(hits, k, input_file, reference):
    row = {}
    row["hit_fraction"] = hits
    row["k"] = k
    row["input_file"] = input_file
    row["reference"] = reference
    return(row)


k_list = range(8, 15)
kmer_dicts = [make_kmer_dictionary(refs, k) for k in k_list]
row_list = []
row_list_poly = []
for p in partis_files:
    mutations = process_partis(os.path.join(data_directory, p))
    mutations_poly = process_partis_poly(os.path.join(data_directory, p),
                                         max_spacing=8, max_mutation_rate=.3)
    imf_out = [indexed_motif_finder(mutations, kmer_dict, k) for
               (kmer_dict, k) in zip(kmer_dicts, k_list)]
    pmf_out = [indexed_motif_finder(mutations_poly, kmer_dict, k) for
               (kmer_dict, k) in zip(kmer_dicts, k_list)]
    
    for (imf_out, pmf_out, k) in zip(imf_out,
                                     pmf_out,
                                     k_list):
        row_list.append(make_row(hit_fraction(imf_out),
                                 k,
                                 p,
                                 "imgt_ighv_human"))
        row_list_poly.append(make_row(hit_fraction(pmf_out),
                                      k,
                                      p,
                                      "imgt_ighv_human"))

df = pd.DataFrame(row_list)
df_poly = pd.DataFrame(row_list_poly)
df.to_csv(output_csv)
df_poly.to_csv(output_csv_poly)
