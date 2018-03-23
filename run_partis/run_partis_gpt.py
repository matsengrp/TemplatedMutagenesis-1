#! /usr/bin/env python
# -*- coding: utf-8 -*-

import os
import os.path
import sys
import time
import numpy as np
import subprocess
import getopt

from os.path import join

# parse input arguments
input_directory = ""
output_directory = ""
options, remainder = getopt.getopt(sys.argv[1:], "i:o:")
for (opt, arg) in options:
    if opt == "-i":
        input_directory = arg
    if opt == "-o":
        output_directory = arg


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

SCRATCH_DIR = 'run_partis/_tmp/'
GERMLINE_GPT = 'run_partis/gl_gpt'
SUFFIX_PATH = 'run_partis/seqs_with_suffix'
#PARTIS = 'partis/bin/partis'
PARTIS = '/Users/juliefukuyama/GitHub/partis/bin/partis'


# add the fake D and J genes
def add_suffix(dataset):
    if not os.path.exists(SUFFIX_PATH):
        print "  creating directories"
        os.makedirs(SUFFIX_PATH)
    cmd = ['seqmagick',
           'convert',
           '--apply-function',
           'run_partis/add_suffix.py:add_suffix',
           os.path.join(input_directory, dataset + '_atleast-2.fastq'),
           os.path.join(SUFFIX_PATH, dataset + '_added_suffix.fasta')]
    print 'calling: ', ' '.join(cmd)
    subprocess.call(map(str, cmd))


# run partis
def annotate(dataset, outdir, gl_dir):
    if not os.path.exists(outdir):
        print "  creating directories"
        os.makedirs(outdir)

    cmd = [PARTIS,
           'annotate',
           '--infname',
           os.path.join(SUFFIX_PATH, dataset + '_added_suffix.fasta'),
           '--outfname',
           join(outdir, dataset + '_annotations.csv'),
           '--initial-germline-dir',
           gl_dir,
           '--only-smith-waterman',
           '--dont-write-parameters',
           '--sw-cachefname',
           SCRATCH_DIR + str(time.time() + np.random.randint(10000)),
           '--gap-open-penalty',
           '60']
    print "  calling:", " ".join(cmd)
    with open(join(SCRATCH_DIR, 'annotation_log.txt'), 'w') as f:
        subprocess.call(map(str, cmd), stdout=f, stderr=f)


for dataset in DATASETS:
    add_suffix(dataset)
    # run partis on against the gpt germline gene set
    annotate(dataset, output_directory, GERMLINE_GPT)
