#!/usr/bin/env python
import subprocess
import pysam
from Bio import SeqIO
import os
from string_compare import get_match_properties
from write_query_fq import write_query
from merge_bowtie import merge_bowtie

seed_length = 4
reference = "test_data/lambda_virus.fa"
query_prefix = "test_data/seed_test_query"
output_prefix = "test_data/seed_test_output"
query_sequence = "ACGCTGGCCATGC"
write_query([query_sequence], ["r1"], [4], [4], query_prefix)
# call the script that runs bowtie
subprocess.call(["./bowtie_exact_match.sh",
                 reference, query_prefix, str(seed_length), output_prefix],
                stdout=open(os.devnull, 'wb'),
                stderr=open(os.devnull, 'wb'))

samfile_fw_right = pysam.AlignmentFile(output_prefix + "_fw_right.bam", "rb")
samfile_fw_left = pysam.AlignmentFile(output_prefix + "_fw_left.bam", "rb")
samfile_rc_right = pysam.AlignmentFile(output_prefix + "_rc_right.bam", "rb")
samfile_rc_left = pysam.AlignmentFile(output_prefix + "_rc_left.bam", "rb")
fw_merged_dict = merge_bowtie(samfile_fw_right, samfile_fw_left, seed_length)
rc_merged_dict = merge_bowtie(samfile_rc_right, samfile_rc_left, seed_length)

# get the match lengths from bowtie
lengths_bt_fw = [len(fw_merged_dict[k]) for k in fw_merged_dict.keys()]
lengths_bt_rc = [len(rc_merged_dict[k]) for k in rc_merged_dict.keys()]

# get the match lengths from the simple function
ref_seq = SeqIO.read(reference, "fasta")
properties_fw = get_match_properties(str(ref_seq.seq), query=query_sequence,
                                     seed_start=4, seed_len=seed_length)
properties_rc = get_match_properties(str(ref_seq.seq), query=query_sequence,
                                     seed_start=4, seed_len=seed_length, rc=True)
lengths_fw = [l for (i, l, left, right) in properties_fw]
lengths_rc = [l for (i, l, left, right) in properties_rc]

print sorted(lengths_bt_fw) == sorted(lengths_fw)
print sorted(lengths_bt_rc) == sorted(lengths_rc)
