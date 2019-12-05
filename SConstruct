#!/usr/bin/env scons
# -*- python -*-

import os
env = Environment(ENV = os.environ)
DATA_DIR = 'data'
OUTPUT_DIR = 'output'
PARTIS = 'partis/bin/partis'
if not os.path.exists(OUTPUT_DIR):
    os.mkdir(OUTPUT_DIR)

## Mutation calling with partis ##

# gpt sequences
partis_gpt = env.Command(
    os.path.join(OUTPUT_DIR, 'partis/partis_output_gpt'),
    [os.path.join(DATA_DIR, 'yeap/presto_output'), PARTIS],
    'python analysis/run_partis_gpt.py --input-directory ${SOURCES[0]} --partis ${SOURCES[1]} --output-directory $TARGET')

# Ebola sequences
partis_ebola = env.Command(
    os.path.join(OUTPUT_DIR, 'partis/partis_output_ebola'),
    [os.path.join(DATA_DIR, 'ebola/ebola_sequences_heavy.fasta'), PARTIS],
    'python analysis/run_partis_ebola.py --input-file ${SOURCES[0]} --partis ${SOURCES[1]} --output-directory $TARGET')

# VB18 sequences
partis_vb18 = env.Command(
    os.path.join(OUTPUT_DIR, 'partis/partis_output_vb18'),
    [os.path.join(DATA_DIR, 'yeap/presto_output'), PARTIS],
    'python analysis/run_partis_vb18.py --input-directory ${SOURCES[0]} --partis ${SOURCES[1]} --output-directory $TARGET')

## Create the mock gpt donor sets ##

# Create the mock gpt donor set based on the human IMGT V genes
gpt_mock_from_human = env.Command(
    os.path.join(OUTPUT_DIR, 'gpt_mock_from_human.fasta'),
    os.path.join(DATA_DIR, 'reference_sets/scripts/imgt_human_v.tree'),
    'python analysis/simulate_gpt_set.py --input-tree $SOURCE --output-fasta $TARGET'
)

# Create the mock gpt donor set based on the mouse IMGT V genes
gpt_mock_from_mouse = env.Command(
    os.path.join(OUTPUT_DIR, 'gpt_mock_from_mouse.fasta'),
    os.path.join(DATA_DIR, 'reference_sets/scripts/imgt_ighv_mouse.tree'),
    'python analysis/simulate_gpt_set.py --input-tree $SOURCE --output-fasta $TARGET'
)

## Compute naive PolyMF rates ##

# Compute the PolyMF rate on the Ebola sequences with the human IMGT V reference set
pmf_rate_ebola_vs_imgt_v = env.Command(
    os.path.join(OUTPUT_DIR, 'pmf_rate_ebola_imgt_v.csv'),
    [partis_ebola,
     os.path.join(DATA_DIR, 'reference_sets/imgt_ighv_human.fasta')],
    'python analysis/compute_mf_rate.py --input-directory ${SOURCES[0]} --reference-fasta ${SOURCES[1]} --output-csv $TARGET'
)

# Compute the PolyMF rate on the VB1-8 sequences with the mouse IMGT V reference set
pmf_rate_vb18_vs_imgt_v = env.Command(
    os.path.join(OUTPUT_DIR, 'pmf_rate_vb18_vs_imgt_v.csv'),
    [partis_vb18,
     os.path.join(DATA_DIR, 'reference_sets/imgt_ighv_mouse.fasta')],
    'python analysis/compute_mf_rate.py --input-directory ${SOURCES[0]} --reference-fasta ${SOURCES[1]} --output-csv $TARGET'
)

# Compute the PolyMF rate on the gpt sequences with the mouse IMGT V reference set
pmf_rate_gpt_vs_imgt_v = env.Command(
    os.path.join(OUTPUT_DIR, 'pmf_rate_gpt_vs_imgt_v.csv'),
    [partis_gpt,
     os.path.join(DATA_DIR, 'reference_sets/imgt_ighv_mouse.fasta')],
    'python analysis/compute_mf_rate.py --input-directory ${SOURCES[0]} --reference-fasta ${SOURCES[1]} --output-csv TARGET'
)

# Compute the PolyMF rate on the gpt sequences with the mock gpt set based on the human V genes
pmf_rate_gpt_vs_human_mock = env.Command(
    os.path.join(OUTPUT_DIR, 'pmf_rate_gpt_vs_human_mock.csv'),
    [partis_gpt, gpt_mock_from_human],
    'python analysis/compute_mf_rate.py --input-directory ${SOURCES[0]} --reference-fasta ${SOURCES[1]} --output-csv TARGET'
)

# Compute the PolyMF rate on the gpt sequences with the mock gpt set based on the mouse V genes
pmf_rate_gpt_vs_mouse_mock = env.Command(
    os.path.join(OUTPUT_DIR, 'pmf_rate_gpt_vs_mouse_mock.csv'),
    [partis_gpt, gpt_mock_from_mouse],
    'python analysis/compute_mf_rate.py --input-directory ${SOURCES[0]} --reference-fasta ${SOURCES[1]} --output-csv TARGET'
)

## Bounds on the rate of templated mutagenesis ##

# Bound on rate of TM in the ebola sequences
env.Command(
    [os.path.join(OUTPUT_DIR, 'ebola_bound.tex'), os.path.join(OUTPUT_DIR, 'ebola_bound.csv')],
    [pmf_rate_ebola, pmf_rate_gpt_vs_human_mock],
    'Rscript analysis/compute_bound.R --naive-rate ${SOURCES[0]} --gpt-fpr ${SOURCES[1]} --output-tex ${TARGETS[0]} --output-csv ${TARGETS[1]} --dale True'
)

# Bound on rate of TM in the VB1-8 sequences
env.Command(
    [os.path.join(OUTPUT_DIR, 'vb18_bound.tex'), os.path.join(OUTPUT_DIR, 'vb18_bound.csv')],
    [pmf_rate_vb18, pmf_rate_gpt_vs_mouse_mock],
    'Rscript analysis/compute_bound.R --naive-rate ${SOURCES[0]} --gpt-fpr ${SOURCES[1]} --output-tex ${TARGETS[0]} --output-csv ${TARGETS[1]} --dale True'
)

## Probabilities of templated mutagenesis given different donor sets ##

## Plots ##


