import pandas as pd
import numpy as np
from process_partis import process_partis


def seed_starts(idx, seed_len, seq_len):
    """Finds starting positions for windows around mutations

    Keyword arguments:
    idx -- Either a single number (for a single mutation) or a tuple or list 
    containing the indices of multiple mutations.
    seed_len -- The length of the window.
    seq_len -- The length of the sequence the mutations occur in.

    Returns: A pair with the lowest and highest indices a window containing 
    the mutations can start at.

    """
    idx_hi = np.max(idx)
    idx_lo = np.min(idx)
    min_start = max(0, idx_hi - seed_len + 1)
    max_start = min(seq_len - seed_len, idx_lo)
    return((min_start, max_start))


def make_kmer_dictionary(references, k):
    """Make a dictionary mapping kmers to sequences they appear in

    Keyword arguments:
    references -- A list of SeqRecord objects (with names and sequences)
    k -- The size of the k-mer

    Returns: A dictionary keyed by k-mers, mapping to sets of
    (name, location) pairs describing where the k-mer occurred.

    """
    d = {}
    for ref in references:
        seq = str(ref.seq)
        seq_len = len(seq)
        for start in range(0, seq_len - k + 1):
            kmer = seq[start:(start + k)]
            if kmer in d.keys():
                d[kmer].add((ref, start))
            else:
                d[kmer] = set([(ref, start)])
    return d


def n_alignments_per_mutation(mutations, kmer_dict, k):
    """ Find the number of unique alignments in the reference set for each 
    mutation

    Keyword arguments:
    mutations -- A data frame created by process_partis
    kmer_dict -- A kmer dictionary created by make_kmer_dictionary
    k -- The k used in the kmer dictionary

    Returns: A data frame with the query name, mutation index, and the
    number of alignments in the reference set explaining that mutation.
    """

    # a data frame describing the matches for each mutation in the references
    imf = indexed_motif_finder(mutations, kmer_dict, k)
    # dictionary that will hold the mutations and how many templates they have
    count_dict = {}
    # Each row in imf describes a mutation. If there is no template
    # for that mutation, the reference_alignment column is
    # np.nan. Otherwise, each row gives the location of one of the
    # templates for that mutation.
    for index, row in imf.iterrows():
        query = row["query_name"]
        query_index = row["query_mutation_index"]
        if np.isnan(row["reference_alignment"]):
            increment = 0
        else:
            increment = 1
        # mutations are described as a pair with the query name and
        # the location of the mutation
        if (query, query_index) in count_dict.keys():
            count_dict[(query, query_index)] += increment
        else:
            count_dict[(query, query_index)] = increment
    # build a DataFrame with query names, indices, and number of alignments
    rows = []
    for (query, query_index) in count_dict.keys():
        rows.append({
                "query_mutation_index": query_index,
                "query_name": query,
                "n_alignments": count_dict[(query, query_index)]
                })
    return(pd.DataFrame(rows))


def indexed_motif_finder(mutations, kmer_dict, k):
    """Find matches around a set of mutations.

    Keyword arguments:
    mutations -- A data frame containing the mutated sequences and
    mutation indices.
    kmer_dict -- A dictionary indexed by k-mers giving the sequences
    they appear in.
    k -- The k used in the kmer dictionary

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
        # make the query sequence by replacing the naive sequence with
        # the mutation we are searching for
        sequence_list = list(row["naive_seq"])
        # for motif_finder/single mutations, we replace just the one base
        if type(row["mutation_index"]) is int:
            sequence_list[row["mutation_index"]] = row["mutated_base"]
        # for poly_motif_finder/multiple mutations, we have to loop
        # over a tuple containing the mutations
        elif type(row["mutation_index"]) is tuple:
            for (i, b) in zip(list(row["mutation_index"]), list(row["mutated_base"])):
                sequence_list[i] = b
        else:
            raise ValueError()
        q = "".join(sequence_list)
        q_id = row["mutated_seq_id"]
        seq_len = len(q)
        mut_idx = row["mutation_index"]
        # search through all the windows around the mutation and check whether 
        # they occur in the references
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
        # if there wasn't a match, we still put the sequence in DataFrame,
        # with np.nan as the value for reference_alignment
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
    """Extends matches from indexed_motif_finder. Adds a columns for left-most 
    and right-most match indices and the match extent"""

    row_count = df.shape[0]
    df["match_extent"] = pd.Series([0] * row_count, index=df.index)
    df["query_left_idx"] = pd.Series([0] * row_count, index=df.index)
    df["query_right_idx"] = pd.Series([0] * row_count, index=df.index)
    for row in range(0, row_count):
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
            elif ref_seq[ref_idx + right + 1] == \
                    query_seq[query_idx + right + 1]:
                right += 1
            else:
                break
        df.loc[row, "match_extent"] = left + right + 1
        df.loc[row, "query_left_idx"] = query_idx - left
        df.loc[row, "query_right_idx"] = query_idx + right


def hit_fraction(df):
    """Fraction of templated mutations.

    Keyword arguments:
    df -- A pandas DataFrame created by indexed_motif_finder.
    Returns: The fraction of mutations in df that have a template.
    """
    # a dictionary, keyed by (query sequence, mutation index) pairs,
    # mapping to 0 if there was no template and 1 if there was a
    # template
    hit_dict = {}
    for (index, row) in df.iterrows():
        query = row["query_sequence"]
        mut_idx = row["query_mutation_index"]
        # if there is no alignment, put a zero for the mutation
        if np.isnan(row["reference_alignment"]):
            hit_dict[(query, mut_idx)] = 0
        # if there is an alignment, put one for the mutation
        else:
            hit_dict[(query, mut_idx)] = 1
    hits = np.mean([hit_dict[k] for k in hit_dict.keys()])
    return hits


def likelihood_given_gcv(partis_file, kmer_dict, k):
    """Finds the likelihood of mutations conditional on being due to gcv

    Keyword arguments:
    partis_file -- A partis csv describing the mutations.
    kmer_dict -- A kmer dictionary describing the references.
    k -- The minimum match length for gcv tracts.

    Returns: A data frame giving the probability of seeing each
    observed mutation. Mutations are described by the name of the
    query sequence and the position of the mutation in that query
    sequence.

    """
    bases = ["A", "C", "G", "T"]
    mut_df = process_partis(partis_file)
    # make a data frame containing all the mutations we didn"t see
    unobs_mut_rows = []
    for index, row in mut_df.iterrows():
        for b in bases:
            if b not in set([row["gl_base"], row["mutated_base"]]):
                r = row.copy()
                unseen_seq = list(r["mutated_seq"])
                unseen_seq[r["mutation_index"]] = b
                r["mutated_seq"] = "".join(unseen_seq)
                r["mutated_base"] = b
                unobs_mut_rows.append(r)
    unobs_mut_df = pd.DataFrame(unobs_mut_rows)

    # run motif finder on the observed and unobserved mutations
    motifs_obs = n_alignments_per_mutation(mut_df, kmer_dict, k)
    motifs_unobs = n_alignments_per_mutation(unobs_mut_df, kmer_dict, k)
    obs_and_unobs = pd.merge(motifs_obs, motifs_unobs,
                             how="outer",
                             on=["query_mutation_index", "query_name"],
                             validate="one_to_one")

    # get the probabilities of seeing the observed
    # mutations. n_alignments_x is the number of alignments for
    # observed mutations because motifs_obs was in the first position
    # in pd.merge
    def get_prob(row):
        n_obs = row["n_alignments_x"]
        n_unobs = row["n_alignments_y"]
        if n_obs + n_unobs == 0:
            return(np.nan)
        return n_obs / (n_obs + n_unobs + 0.)

    obs_and_unobs["prob"] = obs_and_unobs.apply(get_prob, axis=1)
    return(obs_and_unobs)
