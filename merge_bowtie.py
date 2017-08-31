import pysam
from Bio import SeqIO


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
    ref_name = samfile_left.references[0]
    # make a dictionary containing the forward reads indexed by their position
    right_dict = {}
    for read in samfile_right.fetch(ref_name):
        # If the alignment starts with a match of at least seed_length,
        # put it in the dictionary
        if starts_with_match(read.cigartuples, seed_length, left=False):
            right_dict[(read.get_reference_positions()[0], read.query_name)] \
                = read.get_reference_sequence()
    # make a dictionary with the merged reads indexed by their position
    merged_dict = {}
    for read in samfile_left.fetch(ref_name):
        if starts_with_match(read.cigartuples, seed_length, left=True):
            index = read.get_reference_positions()[-seed_length]
            name = read.query_name
            if (index, name) in right_dict.keys():
                right_seq = right_dict[(index, name)]
                merged_seq = read.get_reference_sequence() + right_seq[seed_length:]
                merged_dict[(read.get_reference_positions()[0], name)] = merged_seq
            else:
                # this shouldn't happen
                print("Error: couldn't find a match for left read of query %s at position %i" %
                      (name, read.get_reference_positions()[0]))
                print("read sequence was:")
                print(read.get_reference_sequence())
                print("cigar string was:")
                print(read.cigartuples)
    return(merged_dict)
