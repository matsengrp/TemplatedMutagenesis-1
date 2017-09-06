import pandas as pd


# this should eventually include insertions as well as mutations, for
# now it just compares the naive sequence to the indel reversed
# sequences
def process_partis(partis_file):
    annotations = pd.read_csv(partis_file)
    mutation_map = {}
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
        mutation_map[mutated_seq] = mutation_idx
    return(mutation_map)
