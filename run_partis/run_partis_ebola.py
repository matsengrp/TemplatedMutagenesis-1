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
input_file = ""
output_directory = ""
options, remainder = getopt.getopt(sys.argv[1:], "i:o:")
for (opt, arg) in options:
    if opt == "-i":
        input_file = arg
    if opt == "-o":
        output_directory = arg


SCRATCH_DIR = 'run_partis/_tmp/'
#PARTIS = 'partis/bin/partis'
PARTIS = '/Users/juliefukuyama/GitHub/partis/bin/partis'


# run partis
if not os.path.exists(output_directory):
    print "  creating directories"
    os.makedirs(output_directory)

cmd = [PARTIS,
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
with open(join(SCRATCH_DIR, 'annotation_log.txt'), 'w') as f:
    subprocess.call(map(str, cmd), stdout=f, stderr=f)
