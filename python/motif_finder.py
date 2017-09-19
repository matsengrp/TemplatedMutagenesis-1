from string_compare import extend_one_match
import regex as re
import pandas as pd


# For each (read, mutation) pair, looks at windows of seed_len around
# the mutation and finds all the reference sequences that match.
def motif_finder(mutations, references, seed_len):
    """Find matches around a set of mutations.
    Keyword arguments:
    mutations -- A dictionary with query sequences as keys and a list
    of mutation indices as values.
    references -- A list of reference sequences to match.
    seed_len -- The length of the matched string.
    Returns: A list of Match objects describing the hits.
    """
    match_list = []
    for q in mutations.keys():
        seq_len = len(q)
        for mut_idx in mutations[q]:
            (min_start, max_start) = seed_starts(mut_idx,
                                                 mut_idx,
                                                 seed_len,
                                                 seq_len)
            match_object = Match(q, mut_idx)
            for start in range(min_start, max_start + 1):
                # create the seed
                seed = q[start:(start + seed_len)]
                mut_offset = mut_idx - start
                for r in references:
                    match_indices = [m.start() for m in re.finditer(seed, str(r.seq), overlapped=True)]
                    for match_idx in match_indices:
                        match_object.add_hit(r.seq,
                                             r.name,
                                             match_idx + mut_offset,
                                             seed_len)
            match_list.append(match_object)
    return match_list


def poly_motif_finder(mutations, references, seed_len):
    """Finds matches around pairs of mutations

    Keyword arguments:
    mutations -- A dictionary with query sequences as keys and a list
    of pairs of mutation indices as values.
    references -- A list of the reference sequences to match.
    seed_len -- The length of the matched string.

    Returns: A list of Match objects describing the hits.

    """
    match_list = []
    for q in mutations.keys():
        seq_len = len(q)
        for (idx_lo, idx_hi) in mutations[q]:
            (min_start, max_start) = seed_starts(idx_lo,
                                                 idx_hi,
                                                 seed_len,
                                                 seq_len)
            match_object = Match(q, (idx_lo, idx_hi))
            for start in range(min_start, max_start + 1):
                # create the seed
                seed = q[start:(start + seed_len)]
                offset_lo = idx_lo - start
                offset_hi = idx_hi - start
                for r in references:
                    match_indices = [m.start() for m in re.finditer(seed, str(r.seq), overlapped=True)]
                    for match_idx in match_indices:
                        match_object.add_hit(r.seq,
                                             r.name,
                                             (offset_lo + match_idx,
                                              offset_hi + match_idx),
                                             seed_len)
            match_list.append(match_object)
    return match_list


def seed_starts(idx_lo, idx_hi, seed_len, seq_len):
    """Finds starting positions for windows around mutations

    Keyword arguments:
    idx_lo -- The lowest index of a mutation in the window
    idx_hi -- The highest index of a mutation in the window.
    seed_len -- The length of the window.
    seq_len -- The length of the sequence thu mutations occur in.

    Returns: A pair with the lowest and highest indices the window can
    start at.

    """
    min_start = max(0, idx_hi - seed_len + 1)
    max_start = min(seq_len - seed_len, idx_lo)
    return((min_start, max_start))


def longest_motif_finder(mutations, references, seed_len):

    """For each mutation, finds the longest matches in the references around
    that mutation that match at least min_length bases.

    Keyword arguments:
    mutations -- A dictionary with query sequences as keys and a list
    of mutation indices as values.
    references -- A list of the reference sequences to match.
    seed_len -- The length of the matched string

    Returns: A list of Match objects describing the hits.

    """

    match_list = []
    for q in mutations.keys():
        seq_len = len(q)
        for mut_idx in mutations[q]:
            match_object = Match(q, mut_idx)
            (min_start, max_start) = seed_starts(mut_idx, mut_idx, seed_len, seq_len)
            # get the hits for each window around the mutation
            for r in references:
                for start in range(min_start, max_start + 1):
                    # create the seed
                    seed = q[start:(start + seed_len)]
                    # find all the places the seed matches the reference
                    matches = [m.start() for m in re.finditer(seed, str(r.seq), overlapped=True)]
                    for m in matches:
                        match_extended = extend_one_match(r.seq,
                                                          q,
                                                          start,
                                                          seed_len,
                                                          m)
                        # mut_idx - start + m is the index of the
                        # mutation in the reference sequence, and
                        # match_extended[1] is the length of the
                        # extended match
                        match_object.add_hit(r.seq,
                                             r.name,
                                             mut_idx - start + m,
                                             match_extended[1])
            match_list.append(match_object)
    return match_list


class Match(object):

    def __init__(self, query_seq, mut_idx):
        """Return a new Match object"""
        # the sequence of the query
        self.query = query_seq
        # this can be either a single number or a tuple
        self.mut_idx = mut_idx
        self.ref_seq = []
        self.ref_name = []
        self.ref_idx = []
        self.match_extent = []

    def add_hit(self, ref_seq, ref_name, ref_idx, match_extent):
        self.ref_seq.append(ref_seq)
        self.ref_name.append(ref_name)
        self.ref_idx.append(ref_idx)
        self.match_extent.append(match_extent)

    def num_hits(self):
        return len(self.ref_name)

    def hit_exists(self):
        return len(self.ref_name) > 0


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
                d[kmer].add(ref.name)
            else:
                d[kmer] = set([ref.name])
    return d


def indexed_motif_finder(mutations, kmer_dict, k):
    """Find matches around a set of mutations.

    Keyword arguments:
    mutations -- A dictionary with query sequences as keys and a list
    of mutation indices as values.
    kmer_dict -- A dictionary indexed by k-mers giving the sequences they
    appear in.

    Returns: A pandas DataFrame.

    """
    row_list = []
    for q in mutations.keys():
        seq_len = len(q)
        for mut_idx in mutations[q]:
            (min_start, max_start) = seed_starts(mut_idx,
                                                 mut_idx,
                                                 k,
                                                 seq_len)
            for start in range(min_start, max_start + 1):
                # create the seed
                seed = q[start:(start + k)]
                if seed in kmer_dict:
                    for ref in kmer_dict[seed]:
                        row_list.append({
                            "query_sequence": q,
                            "query_mutation_index": mut_idx,
                            "reference_name": ref
                        })
    return pd.DataFrame(row_list)
