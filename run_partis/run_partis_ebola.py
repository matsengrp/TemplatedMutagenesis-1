#! /usr/bin/env python
# -*- coding: utf-8 -*-

import os
import os.path
import argparse
import subprocess

from os.path import join

# parse input arguments
parser = argparse.ArgumentParser(description='Run partis on ebola sequences')
parser.add_argument('--input-file', dest='input_file')
parser.add_argument('--output-directory', dest='output_directory')
parser.add_argument('--partis', dest='partis')
args = parser.parse_args()


# parse input arguments
input_file = args.input_file
output_directory = args.output_directory
partis = args.partis
scratch_dir = 'run_partis/_tmp/'

# run partis
if not os.path.exists(output_directory):
    print "  creating directories"
    os.makedirs(output_directory)
if not os.path.exists(scratch_dir):
    print "  creating directories"
    os.makedirs(scratch_dir)


cmd = [partis,
       'annotate',
       '--infname',
       input_file,
       '--outfname',
       os.path.join(output_directory, 'ebola_annotations_heavy.csv'),
       '--locus',
       'igh',
       '--species',
       'human'
]
print "  calling:", " ".join(cmd)
with open(join(scratch_dir, 'annotation_log.txt'), 'w') as f:
    subprocess.call(map(str, cmd), stdout=f, stderr=f)
