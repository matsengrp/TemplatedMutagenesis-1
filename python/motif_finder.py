import pandas as pd
import numpy as np


def seed_starts(idx, seed_len, seq_len):
    """Finds starting positions for windows around mutations

    Keyword arguments:
    idx_lo -- The lowest index of a mutation in the window
    idx_hi -- The highest index of a mutation in the window.
    seed_len -- The length of the window.
    seq_len -- The length of the sequence thu mutations occur in.

    Returns: A pair with the lowest and highest indices the window can
    start at.

    """
    idx_hi = np.max(idx)
    idx_lo = np.min(idx)
    min_start = max(0, idx_hi - seed_len + 1)
    max_start = min(seq_len - seed_len, idx_lo)
    return((min_start, max_start))


def make_kmer_dictionary(references, k):
    """Make a dictionary mapping kmers to sequences they appear in

    Keyword arguments:
    references -- A list of SeqRecord objects (have names and sequences)
    k -- The size of the k-mer

    Returns: A dictionary keyed by k-mers, mapping to sets of
    reference names.

    """
    d = {}
    for ref in references:
        seq = str(ref.seq)
        for start in range(0, len(seq) - k + 1):
            kmer = seq[start:(start + k)]
            if kmer in d.keys():
                d[kmer].add((ref, start))
            else:
                d[kmer] = set([(ref, start)])
    return d


def indexed_motif_finder(mutations, kmer_dict, k):
    """Find matches around a set of mutations.

    Keyword arguments:
    mutations -- A dictionary with query sequences as keys and a list
    of mutation indices as values. The list can contain tuples describing
    sets of mutations.
    kmer_dict -- A dictionary indexed by k-mers giving the sequences they
    appear in.

    Returns: A pandas DataFrame containing the query sequence, the
    indices of the mutation(s) in the query, the name and sequence of
    the reference with a match, and a reference alignment. This is the
    position in the reference that matches the mutation in the query
    (if we are loking at single mutations) or the position in the
    reference matching the left-most mutation in a set of mutations in
    the query.

    """
    row_list = []
    for q in mutations.keys():
        seq_len = len(q)
        for mut_idx in mutations[q]:
            (min_start, max_start) = seed_starts(mut_idx,
                                                 k,
                                                 seq_len)
            for start in range(min_start, max_start + 1):
                # create the seed
                seed = q[start:(start + k)]
                mut_offset = np.min(mut_idx) - start
                if seed in kmer_dict:
                    for (ref, ref_idx) in kmer_dict[seed]:
                        row_list.append({
                            "query_sequence": q,
                            "query_mutation_index": mut_idx,
                            "reference_name": ref.name,
                            "reference_sequence": ref.seq,
                            "reference_alignment": ref_idx + mut_offset
                        })
    return pd.DataFrame(row_list)
