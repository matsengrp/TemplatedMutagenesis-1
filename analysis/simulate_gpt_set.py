import pyvolve
import os
import argparse
from Bio import SeqIO

# parse input arguments
parser = argparse.ArgumentParser(description = "Simulate a gpt gene set using pyvolve")
# we need to specify a tree and an output file
parser.add_argument('--input-tree', dest = 'input_tree')
parser.add_argument('--output-fasta', dest = 'output_fasta')
args = parser.parse_args()
fasta_prefix, _ = os.path.splitext(args.output_fasta)
tree_prefix, _ = os.path.splitext(args.input_tree)

# First we need to format the tree for pyvolve --- FastTree's node labels are numbers and pyvolve requires the labels to start with a letter.
# The script also scales up the branch lengths so that the divergences of the simulated sequences match the divergences of the original sequences.
os.system('Rscript analysis/format_tree_for_pyvolve.R --input-tree ' + tree_prefix)
# Read the tree
my_tree = pyvolve.read_tree(file = tree_prefix + '_for_pyvolve.tree')

# Define a codon model, as a pyvolve.Model object, "GY" or "codon" for the GY94-style (uses codon equilibrium frequencies in the matrix)
parameters_alpha_beta = {"beta": 0.65, "alpha": 0.98} # Corresponds to dN/dS = 0.65 / 0.98
my_model = pyvolve.Model("GY", parameters_alpha_beta)
# TODO: fix the last nucleotide in the sequence
my_partition = pyvolve.Partition(models = my_model, root_sequence = "ATGAGCGAAAAATACATCGTCACCTGGGACATGTTGCAGATCCATGCACGTAAACTCGCAAGCCGACTGATGCCTTCTGAACAATGGAAAGGCATTATTGCCGTAAGCCGTGGCGGTCTGGTACCGGGTGCGTTACTGGCGCGTGAACTGGGTATTCGTCATGTCGATACCGTTTGTATTTCCAGCTACGATCACGACAACCAGCGCGAGCTTAAAGTGCTGAAACGCGCAGAAGGCGATGGCGAAGGCTTCATCGTTATTGATGACCTGGTGGATACCGGTGGTACTGCGGTTGCGATTCGTGAAATT")

# Evolve!
my_evolver = pyvolve.Evolver(partitions = my_partition, tree = my_tree)
my_evolver(seqfile = fasta_prefix + '_extra_nucleotide.fasta')

# Cut off the extra nucleotide
f = open(fasta_prefix + '_extra_nucleotide.fasta', 'r')
trimmed_records = []
for record in SeqIO.parse(f, "fasta"):
    trimmed_record = record
    trimmed_record.seq = record.seq[0:308]
    trimmed_records.append(trimmed_record)
## Save the trimmed sequences as fasta
with open(args.output_fasta, 'w') as f:
    SeqIO.write(trimmed_records, f, 'fasta')
