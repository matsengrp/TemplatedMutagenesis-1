from pymotiffinder.motif_finder import likelihood_given_gcv
from pymotiffinder.motif_finder import make_kmer_dictionary,make_kmer_dict_from_fasta
from Bio import SeqIO
import os
import argparse
import pandas as pd
import random

# parse input arguments
parser = argparse.ArgumentParser(description='Compute probability of observed mutations given gcv')
parser.add_argument('--input-directory', dest='input_directory')
parser.add_argument('--output-csv', dest='output_csv')
parser.add_argument('--references', dest='references')
parser.add_argument('--reference-to-drop-kmers', dest='reference_to_drop_kmers', default='')
parser.add_argument('--kmin', dest='kmin', type=int, default=8)
parser.add_argument('--kmax', dest='kmax', type=int, default=14)
parser.add_argument('--max-mutation-rate', dest='max_mutation_rate', type=float, default=.1)
parser.add_argument('--use-indel-seqs', dest='use_indel_seqs', type=bool, default=True)
parser.add_argument('--rc', dest = 'rc', type=bool, default=False)
args = parser.parse_args()

partis_files = os.listdir(args.input_directory)
refs = [r for r in SeqIO.parse(args.references, "fasta")]
reference_name = os.path.splitext(os.path.basename(args.references))[0]
k_list = range(args.kmin, args.kmax + 1)
kmer_dicts = [make_kmer_dictionary(refs, k, reverse_complement=args.rc) for k in k_list]
## drop kmers
if args.reference_to_drop_kmers != "":
    for (i, k) in enumerate(k_list):
        kmer_dict_for_match = make_kmer_dict_from_fasta(args.reference_to_drop_kmers,
                                                        k=k, reverse_complement=args.rc)
        ## The target number of kmers is the size of kmer_dict_for_match
        target_n_kmers = len(kmer_dict_for_match.keys())
        n_kmers_to_remove = len(kmer_dicts[i].keys()) - target_n_kmers
        if n_kmers_to_remove > 0:
            kmers_to_remove = random.sample(kmer_dicts[i].keys(), n_kmers_to_remove)
            for kmer in kmers_to_remove:
                del kmer_dicts[i][kmer]

df_list = []
for f in partis_files:
    for (k, kmer_dict) in zip(k_list, kmer_dicts):
        out = likelihood_given_gcv(os.path.join(args.input_directory, f),
                                   kmer_dict,
                                   k,
                                   max_mutation_rate=args.max_mutation_rate,
                                   use_indel_seqs=args.use_indel_seqs)
        out["k"] = k
        out["source"] = f
        out["reference"] = reference_name
        df_list.append(out)
    print "finished " + f

df = pd.concat(df_list)
df.to_csv(args.output_csv)
