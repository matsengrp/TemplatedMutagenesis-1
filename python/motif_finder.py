

# For each (read, mutation) pair, looks at windows of seed_len around
# the mutation and finds all the reference sequences that match.
def motif_finder(mutations, references, seed_len):
    """Find matches around a set of mutations.

    Keyword arguments:
    mutations -- A dictionary with query sequences as keys and a list
    of mutation indices as values.
    references -- A list of reference sequences to match.
    seed_len -- The length of the matched string.

    Returns: A dictionary keyed by (query, index) tuples with a set of
    reference hits as the values.

    """
    match_dict = {}
    for q in mutations.keys():
        seq_len = len(q)
        for mut_idx in mutations[q]:
            # farthest left starting point for a window of length
            # seed_len containing mut_idx
            min_start = max(0, mut_idx - seed_len + 1)
            # farthest right starting point for a window of length
            # seed_len containing mut_idx
            max_start = min(seq_len - seed_len, mut_idx)
            # a set to store the hits in
            reference_hits = set()
            for start in range(min_start, max_start + 1):
                # create the seed
                seed = q[start:(start + seed_len)]
                for r in references:
                    if seed in r:
                        reference_hits.add(r.name)
            match_dict[(q, mut_idx)] = reference_hits
    return match_dict
