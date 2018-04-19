#!/usr/bin/python

from Bio import SeqIO
from pymotiffinder.motif_finder import make_kmer_dictionary
from pymotiffinder.motif_finder import indexed_motif_finder
from pymotiffinder.motif_finder import hit_fraction
from pymotiffinder.process_partis import process_partis_poly
from pymotiffinder.process_partis import process_partis
import os
import argparse
import pandas as pd

# parse input arguments
parser = argparse.ArgumentParser(description='Run PyMotifFinder')
parser.add_argument('--input-directory', dest='input_directory')
parser.add_argument('--output-csv', dest='output_csv')
parser.add_argument('--output-csv-poly', dest='output_csv_poly')
parser.add_argument('--references', dest='references')
parser.add_argument('--kmin', dest='kmin', type=int, default=8)
parser.add_argument('--kmax', dest='kmax', type=int, default=14)
parser.add_argument('--max-spacing-poly', dest='max_spacing_poly',
                    type=int, default=8)
parser.add_argument('--max-mutation-rate', dest='max_mutation_rate',
                    type=float, default=.3)
args = parser.parse_args()

partis_files = os.listdir(args.input_directory)
refs = [r for r in SeqIO.parse(args.references, "fasta")]
reference_name = os.path.splitext(os.path.basename(args.references))[0]


def make_row(hits, k, input_file, reference):
    row = {}
    row["hit_fraction"] = hits
    row["k"] = k
    row["input_file"] = input_file
    row["reference"] = reference
    return(row)


k_list = range(args.kmin, args.kmax + 1)
kmer_dicts = [make_kmer_dictionary(refs, k) for k in k_list]
row_list = []
row_list_poly = []
for p in partis_files:
    mutations = process_partis(os.path.join(args.input_directory, p))
    mutations_poly = process_partis_poly(
        os.path.join(args.input_directory, p),
        max_spacing=args.max_spacing_poly,
        max_mutation_rate=args.max_mutation_rate)
    imf_out = [indexed_motif_finder(mutations, kmer_dict, k) for
               (kmer_dict, k) in zip(kmer_dicts, k_list)]
    pmf_out = [indexed_motif_finder(mutations_poly, kmer_dict, k) for
               (kmer_dict, k) in zip(kmer_dicts, k_list)]
    
    for (imf_out, pmf_out, k) in zip(imf_out,
                                     pmf_out,
                                     k_list):
        row_list.append(make_row(hit_fraction(imf_out),
                                 k,
                                 os.path.splitext(p)[0],
                                 reference_name))
        row_list_poly.append(make_row(hit_fraction(pmf_out),
                                      k,
                                      os.path.splitext(p)[0],
                                      reference_name))

df = pd.DataFrame(row_list)
df_poly = pd.DataFrame(row_list_poly)
df.to_csv(args.output_csv)
df_poly.to_csv(args.output_csv_poly)
