#!/usr/bin/env python
import subprocess
import pysam
from Bio import SeqIO
import os
from string_compare import get_match_lengths

seed_length = 4
reference = "test_data/seed_test.fa"
query_fw = "test_data/seed_test_query_fw.fq"
query_rc = "test_data/seed_test_query_rc.fq"
# call the script that runs bowtie
subprocess.call(["./bowtie_exact_match.sh",
                 reference, query_fw, query_rc, str(seed_length)],
                stdout=open(os.devnull, 'wb'),
                stderr=open(os.devnull, 'wb'))

ref_seq = SeqIO.read(reference, "fasta")
samfile_fw = pysam.AlignmentFile("output_bt_fw.bam", "rb")
samfile_rc = pysam.AlignmentFile("output_bt_rc.bam", "rb")
ref_name = samfile_fw.references[0]

# make a dictionary containing the forward reads indexed by their position
fw_dict = {}
for fw_read in samfile_fw.fetch(ref_name):
    fw_dict[fw_read.get_reference_positions()[0]] \
        = fw_read.get_reference_sequence()

# make a dictionary with the merged reads indexed by their position
merged_dict = {}
for rc_read in samfile_rc.fetch(ref_name):
    if rc_read.get_reference_positions()[-seed_length] in fw_dict.keys():
        fw_seq = fw_dict[rc_read.get_reference_positions()[-seed_length]]
        merged_seq = rc_read.get_reference_sequence() + fw_seq[seed_length:]
        merged_dict[rc_read.get_reference_positions()[0]] = merged_seq
    else:
        # this should never happen
        print("Error: couldn't find a match for rc read at position %i",
              rc_read.get_reference_positions()[0])

# get the match lengths from bowtie
lengths_bt = [len(merged_dict[k]) for k in merged_dict.keys()]

# get the match lengths from the simple function
lengths = get_match_lengths(str(ref_seq.seq), query="AAAACCCCGGGG",
                            seed_start=4, seed_len=seed_length)

# check that they are the same
if sorted(lengths) == sorted(lengths_bt):
    print("The set of match lengths for the two methods is the same")
else:
    print("The set of match lengths for the two methods is different")

print("Lengths from simple function:")
print(sorted(lengths))
print("Lengths from bowtie2:")
print(sorted(lengths_bt))
