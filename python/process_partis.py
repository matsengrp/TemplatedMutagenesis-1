import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# this should eventually include insertions as well as mutations, for
# now it just compares the naive sequence to the indel reversed
# sequences
def process_partis(partis_file):
    annotations = pd.read_csv(partis_file)
    mutation_rows = []
    row_count = annotations.shape[0]
    for row in range(row_count):
        naive_seq = annotations["naive_seq"][row]
        indel_reversed_seq = annotations["indel_reversed_seqs"][row]
        if pd.isnull(naive_seq):
            continue
        if pd.isnull(indel_reversed_seq):
            mutated_seq = annotations["input_seqs"][row]
        else:
            mutated_seq = indel_reversed_seq
        mutated_seq_len = len(mutated_seq)
        mutation_idx = [i for i in range(mutated_seq_len) if
                        mutated_seq[i] != naive_seq[i]]
        for i in range(mutated_seq_len):
            if(mutated_seq[i] != naive_seq[i]):
                mutation_rows.append({
                        "mutated_seq": mutated_seq,
                        "naive_seq": naive_seq,
                        "mutated_seq_id": annotations["unique_ids"][row],
                        "mutation_index": i,
                        "gl_base": naive_seq[i],
                        "mutated_base": mutated_seq[i]
                        })
    mutation_df = pd.DataFrame(mutation_rows)
    return(mutation_df)


def process_partis_poly(partis_file, min_spacing):
    # first get all the single mutations using process_partis
    mutation_df = process_partis(partis_file)
    poly_mutation_rows = []
    mutation_df_set = set(mutation_df["mutated_seq_id"]) 
    for seq_id in mutation_df_set:
        mutation_df_subset = mutation_df[mutation_df.mutated_seq_id == seq_id]
        poly_mutations = get_pairs(mutation_df_subset["mutation_index"], \
                min_spacing)
        mutated_seq = mutation_df_subset["mutated_seq"][0]
        naive_seq = mutation_df_subset["naive_seq"][0]
        for pm in poly_mutations:
            poly_mutation_rows.append({
                    "mutated_seq": mutated_seq,
                    "naive_seq": naive_seq,
                    "mutated_seq_id": seq_id,
                    "mutation_index": pm,
                    "gl_base": ''.join(naive_seq[i] for i in pm),
                    "mutated_base": ''.join(mutated_seq[i] for i in pm)
                    })
    return(pd.DataFrame(poly_mutation_rows))


def get_pairs(mut_indices, min_spacing):
    pairs = []
    mut_indices_len = len(mut_indices)
    for i in range(0, mut_indices_len - 1):
        for j in range(i+1, mut_indices_len):
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
    row_count = annotations.shape[0]
    for row in range(row_count):
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
                unseen_seq = SeqRecord(Seq("".join(unseen_seq)), \
                        id=annotations["unique_ids"][row])
                mutation_map[unseen_seq] = [i]
    return(mutation_map)
