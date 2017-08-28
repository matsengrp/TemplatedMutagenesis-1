import pysam
from Bio import SeqIO


def merge_bowtie(samfile_right, samfile_left, seed_length):
    ref_name = samfile_left.references[0]
    # make a dictionary containing the forward reads indexed by their position
    right_dict = {}
    for read in samfile_right.fetch(ref_name):
        right_dict[(read.get_reference_positions()[0], read.query_name)] \
            = read.get_reference_sequence()

    # make a dictionary with the merged reads indexed by their position
    merged_dict = {}
    for read in samfile_left.fetch(ref_name):
        index = read.get_reference_positions()[-seed_length]
        name = read.query_name
        if (index, name) in right_dict.keys():
            right_seq = right_dict[(index, name)]
            merged_seq = read.get_reference_sequence() + right_seq[seed_length:]
            merged_dict[(read.get_reference_positions()[0], name)] = merged_seq
        else:
            # this should never happen
            print("Error: couldn't find a match for left read of query %s at position %i" %
                  (name, read.get_reference_positions()[0]))
            print("read sequence was:")
            print(read.get_reference_sequence())
    return(merged_dict)
