import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# this should eventually include insertions as well as mutations, for
# now it just compares the naive sequence to the indel reversed
# sequences
def process_partis(partis_file):
    annotations = pd.read_csv(partis_file)
    mutation_rows = []
    for row in range(annotations.shape[0]):
        naive_seq = annotations["naive_seq"][row]
        indel_reversed_seq = annotations["indel_reversed_seqs"][row]
        if pd.isnull(naive_seq):
            continue
        if pd.isnull(indel_reversed_seq):
            mutated_seq = annotations["input_seqs"][row]
        else:
            mutated_seq = indel_reversed_seq
        mutation_idx = [i for i in range(len(mutated_seq)) if
                        mutated_seq[i] != naive_seq[i]]
        for i in range(len(mutated_seq)):
            if(mutated_seq[i] != naive_seq[i]):
                mutation_rows.append({
                        "mutated_seq": mutated_seq,
                        "mutated_seq_id": annotations["unique_ids"][row],
                        "mutation_index": i,
                        "gl_base": naive_seq[i],
                        "mutated_base": mutated_seq[i]
                        })
    mutation_df = pd.DataFrame(mutation_rows)
    return(mutation_df)


def process_partis_poly(partis_file, min_spacing):
    # first get all the single mutations using process_partis
    mutation_map = process_partis(partis_file)
    # then make a map of all the pairs of mutations within a minimum
    # distance of each other
    poly_mutation_map = {}
    for seq in mutation_map.keys():
        mut_indices = mutation_map[seq]
        poly_mutation_map[seq] = get_pairs(mut_indices, min_spacing)
    return(poly_mutation_map)


def get_pairs(mut_indices, min_spacing):
    pairs = []
    for i in range(0, len(mut_indices) - 1):
        for j in range(i+1, len(mut_indices)):
            if mut_indices[j] - mut_indices[i] <= min_spacing:
                pairs.append((mut_indices[i], mut_indices[j]))
            else:
                continue
    return(pairs)


def unseen_mutations(partis_file):
    annotations = pd.read_csv(partis_file)
    mutation_map = {}
    bases = ['A', 'G', 'C', 'T']
    # loop through the sequences
    for row in range(annotations.shape[0]):
        # if there were indels we use the indel reversed sequence,
        # otherwise we use the input sequence
        naive_seq = annotations["naive_seq"][row]
        indel_reversed_seq = annotations["indel_reversed_seqs"][row]
        if pd.isnull(naive_seq):
            continue
        if pd.isnull(indel_reversed_seq):
            mutated_seq = annotations["input_seqs"][row]
        else:
            mutated_seq = indel_reversed_seq
        # a list of mutation indices for the sequence
        mutation_idx = [i for i in range(len(mutated_seq)) if
                        mutated_seq[i] != naive_seq[i]]
        # for each mutation and each base that wasn't either germline
        # or observed, make a sequence with the unseen base at that
        # index
        for i in mutation_idx:
            for b in bases:
                if b in [mutated_seq[i], naive_seq[i]] :
                    continue
                unseen_seq = list(mutated_seq)
                unseen_seq[i] = b
                unseen_seq = SeqRecord(Seq("".join(unseen_seq)), id=annotations["unique_ids"][row])
                mutation_map[unseen_seq] = [i]
    return(mutation_map)
