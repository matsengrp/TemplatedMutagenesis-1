from motif_finder import n_alignments_per_mutation
import pandas as pd
from process_partis import process_partis

def templates_per_base(partis_file, kmer_dict, k):
    """

    Keyword arguments:
    partis_file -- A partis csv describing the mutations.
    kmer_dict -- A kmer dictionary describing the references.
    k -- The minimum match length for gcv tracts.

    Returns: 

    """
    bases = ["A", "C", "G", "T"]
    mut_df = process_partis(partis_file)
    output = []
    # make a data frame containing all the mutations we didn"t see
    for index, row in mut_df.iterrows():
        for b in bases:
            r = row.copy()
            r["mutated_seq"] = make_mutated_sequence(r["naive_seq"], r["mutation_index"], b)
            n = n_alignments_per_mutation(r, kmer_dict, k)
            output_row = {"gl_base": r["gl_base"],
                          "template_base": b,
                          "mutated_seq_id": r["mutated_seq_id"],
                          "mutation_index": r["mutation_index"],
                          "true_mutation": r["mutated_base"]}
            output.append(output_row)
    return(pd.DataFrame(output))

def make_mutated_sequence(naive, index, base):
    seq_as_list = list(naive)
    seq_as_list[index] = base
    return("".join(seq_as_list))
