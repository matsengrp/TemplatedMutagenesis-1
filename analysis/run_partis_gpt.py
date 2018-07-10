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
parser = argparse.ArgumentParser(description='Run partis on gpt sequences')
parser.add_argument('--input-directory', dest='input_directory')
parser.add_argument('--output-directory', dest='output_directory')
parser.add_argument('--partis', dest='partis')
args = parser.parse_args()


# parse input arguments
input_directory = args.input_directory
output_directory = args.output_directory
partis = args.partis

# the gpt passenger datasets
DATASETS = [
    'SRR2229714',
    'SRR2229715',
    'SRR2229716',
    'SRR2229717',
    'SRR2229718',
    'SRR2229719',
    'SRR2229702',
    'SRR2229703',
    'SRR2229704',
    'SRR2229705',
    'SRR2229706',
    'SRR2229707'
    ]

SCRATCH_DIR = '_tmp/'
GERMLINE_GPT = 'data/reference_sets/gl_gpt'
SUFFIX_PATH = 'seqs_with_suffix'


# add the fake D and J genes
def add_suffix(dataset):
    if not os.path.exists(SUFFIX_PATH):
        print "  creating directories"
        os.makedirs(SUFFIX_PATH)
    cmd = ['seqmagick',
           'convert',
           '--apply-function',
           'analysis/add_suffix.py:add_suffix',
           os.path.join(input_directory, dataset + '_atleast-2.fastq'),
           os.path.join(SUFFIX_PATH, dataset + '_added_suffix.fasta')]
    print 'calling: ', ' '.join(cmd)
    subprocess.call(map(str, cmd))


# run partis
def annotate(dataset, outdir, gl_dir):
    if not os.path.exists(outdir):
        print "  creating directories"
        os.makedirs(outdir)
    if not os.path.exists(SCRATCH_DIR):
        print "  creating directories"
        os.makedirs(SCRATCH_DIR)
    cmd = [partis,
           'annotate',
           '--infname',
           os.path.join(SUFFIX_PATH, dataset + '_added_suffix.fasta'),
           '--outfname',
           join(outdir, dataset + '_annotations.csv'),
           '--initial-germline-dir',
           gl_dir,
           '--only-smith-waterman',
           '--sw-cachefname',
           os.path.join(SCRATCH_DIR, str(time.time() + np.random.randint(10000))),
           '--gap-open-penalty',
           '60']
    print "  calling:", " ".join(cmd)
    with open(join(SCRATCH_DIR, 'annotation_log_' + dataset + '.txt'), 'w') as f:
        subprocess.call(map(str, cmd), stdout=f, stderr=f)


for dataset in DATASETS:
    add_suffix(dataset)
    # run partis on against the gpt germline gene set
    annotate(dataset, output_directory, GERMLINE_GPT)
