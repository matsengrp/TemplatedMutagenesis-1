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
