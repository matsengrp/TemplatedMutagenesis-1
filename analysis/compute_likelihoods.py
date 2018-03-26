from gcgcgc.likelihood_given_gcv import likelihood_given_gcv
from gcgcgc.motif_finder import make_kmer_dictionary
from Bio import SeqIO
import os
import sys
import getopt
import pandas as pd

data_directory = ""
output_csv = ""
options, remainder = getopt.getopt(sys.argv[1:], "i:o:")
for (opt, arg) in options:
    if opt == "-i":
        data_directory = arg
    if opt == "-o":
        output_csv = arg

reference_file_gpt = "./data/reference_sets/gpt_132.fasta"
reference_file_mus = "./data/reference_sets/mus_musculus_129S1_v_genes.fasta"

partis_files = os.listdir(data_directory)
k_list = range(8, 15)
references_gpt = [r for r in SeqIO.parse(reference_file_gpt, "fasta")]
references_mus = [r for r in SeqIO.parse(reference_file_mus, "fasta")]
kmer_dicts_gpt = [make_kmer_dictionary(references_gpt, k) for k in k_list]
kmer_dicts_mus = [make_kmer_dictionary(references_mus, k) for k in k_list]
df_list = []
for f in partis_files:
    for (k, kmer_dict_gpt, kmer_dict_mus) in zip(k_list, kmer_dicts_gpt, kmer_dicts_mus):
        out_gpt = likelihood_given_gcv(os.path.join(data_directory, f), kmer_dict_gpt, k)
        out_mus = likelihood_given_gcv(os.path.join(data_directory, f), kmer_dict_mus, k)
        out_gpt["k"] = k
        out_mus["k"] = k
        out_gpt["source"] = f
        out_mus["source"] = f
        out_gpt["reference"] = "gpt"
        out_mus["reference"] = "mus"
        df_list.append(out_gpt)
        df_list.append(out_mus)
    print "finished " + f

df = pd.concat(df_list)
df.to_csv(output_csv)
