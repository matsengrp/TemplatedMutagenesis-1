from motif_finder import likelihood_given_gcv
from motif_finder import make_kmer_dictionary
from Bio import SeqIO
import os
import argparse
import pandas as pd

# parse input arguments
parser = argparse.ArgumentParser(description='Compute probability of observed mutations given gcv')
parser.add_argument('--input-directory', dest='input_directory')
parser.add_argument('--output-csv', dest='output_csv')
parser.add_argument('--references', dest='references')
parser.add_argument('--kmin', dest='kmin', type=int, default=8)
parser.add_argument('--kmax', dest='kmax', type=int, default=14)
args = parser.parse_args()

partis_files = os.listdir(args.input_directory)
refs = [r for r in SeqIO.parse(args.references, "fasta")]
reference_name = os.path.splitext(os.path.basename(args.references))[0]
k_list = range(args.kmin, args.kmax)
kmer_dicts = [make_kmer_dictionary(refs, k) for k in k_list]
df_list = []
for f in partis_files:
    for (k, kmer_dict) in zip(k_list, kmer_dicts):
        out = likelihood_given_gcv(os.path.join(args.input_directory, f), kmer_dict, k)
        out["k"] = k
        out["source"] = f
        out["reference"] = reference_name
        df_list.append(out)
    print "finished " + f

df = pd.concat(df_list)
df.to_csv(args.output_csv)
