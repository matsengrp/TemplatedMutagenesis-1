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


def n_alignments_per_mutation(mutations, kmer_dict, k):
    """ Find the number of unique alignments in the reference set for each mutation

    Keyword arguments:
    mutations -- A data frame created by process_partis
    kmer_dict -- A kmer dictionary created by make_kmer_dictionary
    k -- The k used in the kmer dictionary

    Returns: A data frame with the query name, mutation index, and the
    number alignments in the reference set explaining that mutation.
    """

    # find all the matches in the references
    imf = indexed_motif_finder(mutations, kmer_dict, k)
    # extends the matches and drops all the duplicate matches from
    # overlapping windows
    extend_matches(imf)
    # data frame with the unique mutations
    query_and_idx = imf[["query_mutation_index", "query_name"]].drop_duplicates()
    rows = []
    for (q, i) in zip(query_and_idx["query_name"], query_and_idx["query_mutation_index"]):
        templates = imf.loc[(imf.query_name == q) & (imf.query_mutation_index == i), :]
        # no alignment is encoded by nan in reference_alignment
        n_align = len([a for a in templates["reference_alignment"] if not np.isnan(a)])
        rows.append({
                "query_mutation_index": i,
                "query_name": q,
                "n_alignments": n_align
                })
    return(pd.DataFrame(rows))


def indexed_motif_finder(mutations, kmer_dict, k):
    """Find matches around a set of mutations.

    Keyword arguments: 
    mutations -- A data frame containing the mutated sequences and
    mutation indices.
    kmer_dict -- A dictionary indexed by k-mers giving the sequences
    they appear in.

    Returns: A pandas DataFrame containing the query sequence, the
    indices of the mutation(s) in the query, the name and sequence of
    the reference with a match, and a reference alignment. This is the
    position in the reference that matches the mutation in the query
    (if we are loking at single mutations) or the position in the
    reference matching the left-most mutation in a set of mutations in
    the query.

    """
    row_list = []
    for index, row in mutations.iterrows():
        q = row["mutated_seq"]
        q_id = row["mutated_seq_id"]
        seq_len = len(q)
        mut_idx = row["mutation_index"]
        found_match = False
        (min_start, max_start) = seed_starts(mut_idx, k, seq_len)
        for start in range(min_start, max_start + 1):
            # create the seed
            seed = q[start:(start + k)]
            mut_offset = np.min(mut_idx) - start
            if seed in kmer_dict:
                for (ref, ref_idx) in kmer_dict[seed]:
                    row_list.append({
                            "query_sequence": q,
                            "query_name": q_id,
                            "query_mutation_index": mut_idx,
                            "reference_name": ref.name,
                            "reference_sequence": str(ref.seq),
                            "reference_alignment": ref_idx + mut_offset
                            })
                found_match = True
        if not found_match:
            row_list.append({
                    "query_sequence": q,
                    "query_name": q_id,
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


def hit_fraction(df):
    hit_dict = {}
    for row in df.index:
        query = df.loc[row, "query_sequence"]
        mut_idx = df.loc[row, "query_mutation_index"]
        if df.loc[row, "reference_name"] == "":
            hit_dict[(query, mut_idx)] = 0
        else:
            hit_dict[(query, mut_idx)] = 1
    hits = np.mean([hit_dict[k] for k in hit_dict.keys()])
    return hits
