from pymotiffinder.motif_finder import motif_finder
from pymotiffinder.motif_finder import poly_motif_finder
from pymotiffinder.motif_finder import hit_fraction
from pymotiffinder.motif_finder import make_kmer_dict_from_fasta
import os
import argparse
import pandas as pd

# parse input arguments
parser = argparse.ArgumentParser(description='Run PyMotifFinder')
parser.add_argument('--input-directory', dest='input_directory')
parser.add_argument('--output-csv', dest='output_csv')
parser.add_argument('--reference-fasta', dest='reference_fasta')
parser.add_argument('--kmin', dest='kmin', type=int, default=8)
parser.add_argument('--kmax', dest='kmax', type=int, default=14)
parser.add_argument('--max-spacing-poly', dest='max_spacing_poly',
                    type=int, default=8)
parser.add_argument('--max-mutation-rate', dest='max_mutation_rate',
                    type=float, default=.3)
args = parser.parse_args()
partis_files = os.listdir(args.input_directory)
reference_name = os.path.splitext(os.path.basename(args.reference_fasta))[0]


def make_row(hits, k, input_file, reference, type, reverse_complement):
    row = {}
    row["hit_fraction"] = hits
    row["k"] = k
    row["input_file"] = input_file
    row["reference"] = reference
    row["type"] = type
    row["reverse_complement"] = reverse_complement
    return(row)


k_list = range(args.kmin, args.kmax + 1)
row_list = []
row_list_poly = []
for rc in [True, False]:
    for k in k_list:
        kmer_dict = make_kmer_dict_from_fasta(args.reference_fasta,
                                              k=k, reverse_complement=rc)
        for p in partis_files:
            mf_out = motif_finder(
                partis_file=os.path.join(args.input_directory, p),
                reference_fasta=None,
                k=k,
                reverse_complement=rc,
                max_mutation_rate=args.max_mutation_rate,
                kmer_dict=kmer_dict)
            pmf_out = poly_motif_finder(
                partis_file=os.path.join(args.input_directory, p),
                reference_fasta=None,
                k=k,
                max_spacing=8,
                reverse_complement=rc,
                max_mutation_rate=args.max_mutation_rate,
                kmer_dict=kmer_dict)

            row_list.append(make_row(hit_fraction(mf_out),
                                     k,
                                     os.path.splitext(p)[0],
                                     reference_name,
                                     "mf",
                                     rc))
            row_list.append(make_row(hit_fraction(pmf_out),
                                     k,
                                     os.path.splitext(p)[0],
                                     reference_name,
                                     "pmf",
                                     rc))

df = pd.DataFrame(row_list)
df.to_csv(args.output_csv)
