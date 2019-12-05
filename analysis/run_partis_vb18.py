#! /usr/bin/env python
# -*- coding: utf-8 -*-

import os
import os.path
import time
import numpy as np
import subprocess
import argparse

from os.path import join

# parse input arguments
parser = argparse.ArgumentParser(description='Run partis on vb18 sequences')
parser.add_argument('--input-directory', dest='input_directory')
parser.add_argument('--output-directory', dest='output_directory')
parser.add_argument('--partis', dest='partis')
args = parser.parse_args()


# parse input arguments
input_directory = args.input_directory
output_directory = args.output_directory
partis = args.partis

# Where the dataset numbers come from:
# t is the Sra Run table, data/yeap/SraRunTable.txt
# subset(t, genotype_s == "VB18 passenger mice" & allele_s == "Passenger" & tissue_s %in% c("PP GC B cells", "Splenic GC B cells") & immunization_s == "NP-CGG")

# the VB18 passenger datasets
DATASETS = [
    'SRR2229682',
    'SRR2229683',
    'SRR2229684',
    'SRR2229685',
    'SRR2229686',
    'SRR2229687',
    'SRR2229670',
    'SRR2229671',
    'SRR2229672',
    'SRR2229673',
    'SRR2229674',
    'SRR2229675'
]


SCRATCH_DIR = '_tmp/'
FASTA_PATH = 'seqs_as_fasta'

def convert_to_fasta(dataset):
    if not os.path.exists(FASTA_PATH):
        print "  creating directories"
        os.makedirs(FASTA_PATH)
    cmd = ['seqmagick',
           'convert',
           os.path.join(input_directory, dataset + '_atleast-2.fastq'),
           os.path.join(FASTA_PATH, dataset + '_atleast-2.fasta')]
    print 'calling: ', ' '.join(cmd)
    subprocess.call(cmd)


# run partis
def annotate(dataset, outdir):
    if not os.path.exists(outdir):
        print "  creating directories"
        os.makedirs(outdir)
    if not os.path.exists(SCRATCH_DIR):
        print "  creating directories"
        os.makedirs(SCRATCH_DIR)
    cmd = [partis,
           'annotate',
           '--infname',
           os.path.join(FASTA_PATH, dataset + '_atleast-2.fasta'),
           '--outfname',
           os.path.join(outdir, dataset + '_annotations.csv'),
           '--locus',
           'igh',
           '--species',
           'mouse']
    print "  calling:", " ".join(cmd)
    with open(join(SCRATCH_DIR, 'annotation_log_' + dataset + '.txt'), 'w') as f:
        subprocess.call(cmd, stdout=f, stderr=f)


for dataset in DATASETS:
    # run partis against the mouse germline gene set
    convert_to_fasta(dataset)
    annotate(dataset, output_directory)
