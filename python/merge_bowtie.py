import pysam


def starts_with_match(cigartuples, seed_length, left):
    if left:
        tuple_0 = cigartuples[-1]
    else:
        tuple_0 = cigartuples[0]
    if tuple_0[0] == 0 and tuple_0[1] >= seed_length:
        return True
    else:
        return False


def merge_bowtie(samfile_right, samfile_left, seed_length):
    # make a dictionary containing the forward reads indexed by their position
    right_dict = {}
    for read in samfile_right.fetch():
        # If the alignment starts with a match of at least seed_length,
        # put it in the dictionary
        if starts_with_match(read.cigartuples, seed_length, left=False):
            right_dict[(read.reference_name, read.get_reference_positions()[0], read.query_name)] \
                = read.get_reference_sequence()
    # make a dictionary with the merged reads indexed by their position
    merged_dict = {}
    for read in samfile_left.fetch():
        if starts_with_match(read.cigartuples, seed_length, left=True):
            index = read.get_reference_positions()[-seed_length]
            name = read.query_name
            ref_name = read.reference_name
            if (ref_name, index, name) in right_dict.keys():
                right_seq = right_dict[(ref_name, index, name)]
                merged_seq = read.get_reference_sequence() + right_seq[seed_length:]
                merged_dict[(ref_name, read.get_reference_positions()[0], name)] = merged_seq
    return(merged_dict)
