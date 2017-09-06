

# For each (read, mutation) pair, looks at windows of seed_len around
# the mutation and finds all the reference sequences that match.
def sliding_window_match(mutations, references, seed_len):
    match_dict = {}
    for q in mutations.keys():
        seq_len = len(q)
        for mut_idx in mutations[q]:
            min_start = max(0, mut_idx - seed_len + 1)
            max_start = min(seq_len - seed_len + 1, mut_idx)
            seed_to_refs = {}
            for start in range(min_start, max_start + 1):
                seed = q[start:(start + seed_len)]
                seed_to_refs[seed] = get_matches(seed, references)
            match_dict[(q, mut_idx)] = seed_to_refs
    return match_dict


def get_matches(seed, references):
    matches = []
    for r in references:
        if seed in r:
            matches.append(r.name)
        if seed in r.reverse_complement():
            matches.append(r.name)
    return(matches)
