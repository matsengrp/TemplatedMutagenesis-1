import regex as re
from Bio.Seq import Seq


def get_match_properties(reference, query, seed_start, seed_len, rc=False):
    """Finds matches between a seed and a reference sequence, extends the
    seed in both directions, and returns the length of the extended
    match.
    
    Keyword arguments:
    reference -- the reference sequence
    query -- the query sequence
    seed_start -- the position of the seed in the query sequence
    seed_len -- the length of the seed
    rc -- should the matching be on the reverse complement?

    """
    if rc:
        query = query[::-1]
        reference = str(Seq(reference).complement())
        seed_start = len(query) - seed_start - seed_len
    matches = [m.start() for m in re.finditer(query[seed_start:(seed_start + seed_len)], reference, overlapped=True)]
    match_properties = [extend_one_match(reference, query, seed_start, seed_len, m) for m in matches]
    return match_properties


def extend_one_match(reference, query, seed_start, seed_len, match_idx):
    """Extend a match between a seed in the query and a reference
    sequence.
    
    Keyword arguments: 
    reference -- the reference sequence
    query -- the query sequence
    seed_start -- the position of the seed in the query sequence
    seed_len -- the length of the seed in the query sequences
    match_idx -- a position in the reference sequence where the seed matches

    """
    left = 0
    right = 0
    # walk left along the reference and query until the first mismatch
    # or we reach the end
    while True:
        if(match_idx - left - 1 < 0):
            break
        elif(seed_start - left - 1 < 0):
            break
        elif(reference[match_idx - left - 1] == query[seed_start - left - 1]):
            left += 1
        else:
            break
    # walk right along the reference and query until the first
    # mismatch or we reach the end
    while True:
        if(match_idx + seed_len + right >= len(reference)):
            break
        elif(seed_start + seed_len + right >= len(query)):
            break
        elif(reference[match_idx + seed_len + right] == query[seed_start + seed_len + right]):
            right += 1
        else:
            break
    return (match_idx - left, left + right + seed_len, left, right)
