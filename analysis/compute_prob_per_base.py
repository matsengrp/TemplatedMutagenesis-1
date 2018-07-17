from pymotiffinder.motif_finder import make_kmer_dictionary
from pymotiffinder.motif_finder import per_base_alignments
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
parser.add_argument('--max-mutation-rate', dest='max_mutation_rate', type=float, default=.1)
parser.add_argument('--use-indel-seqs', dest='use_indel_seqs', type=bool, default=True)
args = parser.parse_args()

partis_files = os.listdir(args.input_directory)
refs = [r for r in SeqIO.parse(args.references, "fasta")]
reference_name = os.path.splitext(os.path.basename(args.references))[0]
k_list = range(args.kmin, args.kmax + 1)
kmer_dicts = [make_kmer_dictionary(refs, k) for k in k_list]
k_list = range(args.kmin, args.kmax + 1)
refs = [r for r in SeqIO.parse(args.references, "fasta")]
kmer_dicts = [make_kmer_dictionary(refs, k) for k in k_list]
df_list = []
for f in partis_files:
    for (k, kmer_dict) in zip(k_list, kmer_dicts):
        out = per_base_alignments(os.path.join(args.input_directory, f),
                                  kmer_dict, k,
                                  max_mutation_rate=args.max_mutation_rate,
                                  use_indel_seqs=args.use_indel_seqs)
        out["k"] = k
        out["source"] = f
        out["reference"] = reference_name
        df_list.append(out)
    print "finished " + f

df = pd.concat(df_list)
df.to_csv(args.output_csv)
