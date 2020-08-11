from pymotiffinder.motif_finder import poly_motif_finder
from pymotiffinder.motif_finder import make_kmer_dict_from_fasta
import os
import argparse
import pandas as pd
import random

# parse input arguments
parser = argparse.ArgumentParser(description='Run PyMotifFinder')
parser.add_argument('--input-directory', dest='input_directory')
parser.add_argument('--output-csv', dest='output_csv')
parser.add_argument('--reference-fasta', dest='reference_fasta')
parser.add_argument('--reference-to-drop-kmers', dest='reference_to_drop_kmers', default = "")
parser.add_argument('--kmin', dest='kmin', type=int, default=8)
parser.add_argument('--kmax', dest='kmax', type=int, default=14)
parser.add_argument('--max-mutation-rate', dest='max_mutation_rate',
                    type=float, default=.3)
args = parser.parse_args()
partis_files = os.listdir(args.input_directory)
reference_name = os.path.splitext(os.path.basename(args.reference_fasta))[0]


def make_row(hits, k, input_file, reference, type, reverse_complement, dale_method):
    row = {}
    row["hit_fraction"] = hits
    row["k"] = k
    row["input_file"] = input_file
    row["reference"] = reference
    row["type"] = type
    row["reverse_complement"] = reverse_complement
    row["dale_method"] = dale_method
    return(row)

k_list = range(args.kmin, args.kmax + 1)
row_list = []
for rc in [False, True]:
    for dale_method in [True]:
        for k in k_list:
            ## kmer dictionary for donor set
            kmer_dict = make_kmer_dict_from_fasta(args.reference_fasta,
                                                  k=k, reverse_complement=rc)
            ## kmer dictionary for the set we want to match.
            if(args.reference_to_drop_kmers != ""):
                kmer_dict_for_match = make_kmer_dict_from_fasta(args.reference_to_drop_kmers,
                                                                k=k, reverse_complement=rc)
                ## The target number of kmers is the size of kmer_dict_for_match
                target_n_kmers = len(kmer_dict_for_match.keys())
                n_kmers_to_remove = len(kmer_dict.keys()) - target_n_kmers
                if n_kmers_to_remove > 0:
                    kmers_to_remove = random.sample(kmer_dict.keys(), n_kmers_to_remove)
                    for kmer in kmers_to_remove:
                        del kmer_dict[kmer]
            for p in partis_files:
                pmf_out, pmf_frac = poly_motif_finder(
                    partis_file=os.path.join(args.input_directory, p),
                    reference_fasta=None,
                    k=k,
                    reverse_complement=rc,
                    max_mutation_rate=args.max_mutation_rate,
                    kmer_dict=kmer_dict, dale_method = dale_method)
                row_list.append(make_row(pmf_frac,
                                     k,
                                     os.path.splitext(p)[0],
                                     reference_name,
                                     "pmf",
                                     rc, dale_method))

df = pd.DataFrame(row_list)
df.to_csv(args.output_csv)
