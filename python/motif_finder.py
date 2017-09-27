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
            found_match = False
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
                            "reference_sequence": str(ref.seq),
                            "reference_alignment": ref_idx + mut_offset
                        })
                    found_match = True
            if not found_match:
                row_list.append({
                    "query_sequence": q,
                    "query_mutation_index": mut_idx,
                    "reference_name": "",
                    "reference_sequence": "",
                    "reference_alignment": np.nan
                })
    return pd.DataFrame(row_list).drop_duplicates().reset_index(drop=True)


def extend_matches(df):
    """Extends matches from indexed_motif_finder."""

    df["match_extent"] = pd.Series([0] * df.shape[0], index=df.index)
    df["query_left_idx"] = pd.Series([0] * df.shape[0], index=df.index)
    df["query_right_idx"] = pd.Series([0] * df.shape[0], index=df.index)
    for row in range(0, df.shape[0]):
        if df.loc[row, "reference_sequence"] == "":
            continue
        left = 0
        right = 0
        query_idx = int(np.min(df.loc[row, "query_mutation_index"]))
        ref_idx = int(df.loc[row, "reference_alignment"])
        query_seq = df.loc[row, "query_sequence"]
        ref_seq = df.loc[row, "reference_sequence"]
        while True:
            if query_idx - left - 1 < 0:
                break
            elif ref_idx - left - 1 < 0:
                break
            elif ref_seq[ref_idx - left - 1] == query_seq[query_idx - left - 1]:
                left += 1
            else:
                break
        while True:
            if ref_idx + right + 1 >= len(ref_seq):
                break
            elif query_idx + right + 1 >= len(query_seq):
                break
            elif ref_seq[ref_idx + right + 1] == query_seq[query_idx + right + 1]:
                right += 1
            else:
                break
        df.loc[row, "match_extent"] = left + right + 1
        df.loc[row, "query_left_idx"] = query_idx - left
        df.loc[row, "query_right_idx"] = query_idx + right
