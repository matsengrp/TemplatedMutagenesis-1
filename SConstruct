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

## Compute naive PolyMF rates ##

# Compute the PolyMF rate on the Ebola sequences with the human IMGT V reference set
pmf_rate_ebola_vs_human_imgt_v = env.Command(
    os.path.join(OUTPUT_DIR, 'pmf_rate_ebola_human_imgt_v.csv'),
    [partis_ebola,
     os.path.join(DATA_DIR, 'reference_sets/imgt_ighv_human.fasta')],
    'python analysis/compute_mf_rate.py --input-directory ${SOURCES[0]} --reference-fasta ${SOURCES[1]} --output-csv $TARGET'
)

# Compute the PolyMF rate on the VB1-8 sequences with the mouse IMGT V reference set
pmf_rate_vb18_vs_mouse_imgt_v = env.Command(
    os.path.join(OUTPUT_DIR, 'pmf_rate_vb18_vs_mouse_imgt_v.csv'),
    [partis_vb18,
     os.path.join(DATA_DIR, 'reference_sets/imgt_ighv_mouse.fasta')],
    'python analysis/compute_mf_rate.py --input-directory ${SOURCES[0]} --reference-fasta ${SOURCES[1]} --output-csv $TARGET'
)

# Compute the PolyMF rate on the gpt sequences with the mouse IMGT V reference set
pmf_rate_gpt_vs_mouse_imgt_v = env.Command(
    os.path.join(OUTPUT_DIR, 'pmf_rate_gpt_vs_mouse_imgt_v.csv'),
    [partis_gpt,
     os.path.join(DATA_DIR, 'reference_sets/imgt_ighv_mouse.fasta')],
    'python analysis/compute_mf_rate.py --input-directory ${SOURCES[0]} --reference-fasta ${SOURCES[1]} --output-csv $TARGET'
)

# Compute the PolyMF rate on the gpt sequences with the mock gpt set based on the human V genes
pmf_rate_gpt_vs_human_mock = env.Command(
    os.path.join(OUTPUT_DIR, 'pmf_rate_gpt_vs_human_mock.csv'),
    [partis_gpt,
     os.path.join(DATA_DIR, 'reference_sets/gpt_mock_from_human.fasta')],
    'python analysis/compute_mf_rate.py --input-directory ${SOURCES[0]} --reference-fasta ${SOURCES[1]} --output-csv $TARGET'
)

# Compute the PolyMF rate on the gpt sequences with the mock gpt set based on the mouse V genes
pmf_rate_gpt_vs_mouse_mock = env.Command(
    os.path.join(OUTPUT_DIR, 'pmf_rate_gpt_vs_mouse_mock.csv'),
    [partis_gpt,
     os.path.join(DATA_DIR, 'reference_sets/gpt_mock_from_mouse.fasta')],
    'python analysis/compute_mf_rate.py --input-directory ${SOURCES[0]} --reference-fasta ${SOURCES[1]} --output-csv $TARGET'
)

## Bounds on the rate of templated mutagenesis ##

# Bound on rate of TM in the ebola sequences
env.Command(
    [os.path.join(OUTPUT_DIR, 'ebola_bound.tex'), os.path.join(OUTPUT_DIR, 'ebola_bound.csv')],
    [pmf_rate_ebola_vs_human_imgt_v, pmf_rate_gpt_vs_human_mock],
    'Rscript analysis/compute_bound.R --naive-rate ${SOURCES[0]} --gpt-fpr ${SOURCES[1]} --output-tex ${TARGETS[0]} --output-csv ${TARGETS[1]} --dale True'
)

# Bound on rate of TM in the VB1-8 sequences
env.Command(
    [os.path.join(OUTPUT_DIR, 'vb18_bound.tex'), os.path.join(OUTPUT_DIR, 'vb18_bound.csv')],
    [pmf_rate_vb18_vs_mouse_imgt_v, pmf_rate_gpt_vs_mouse_mock],
    'Rscript analysis/compute_bound.R --naive-rate ${SOURCES[0]} --gpt-fpr ${SOURCES[1]} --output-tex ${TARGETS[0]} --output-csv ${TARGETS[1]} --dale True'
)

# Bound on rate of TM in the ebola sequences with reverse complements included in donor set
env.Command(
    [os.path.join(OUTPUT_DIR, 'ebola_bound_rc.tex'), os.path.join(OUTPUT_DIR, 'ebola_bound_rc.csv')],
    [pmf_rate_ebola_vs_human_imgt_v, pmf_rate_gpt_vs_human_mock],
    'Rscript analysis/compute_bound.R --naive-rate ${SOURCES[0]} --gpt-fpr ${SOURCES[1]} --output-tex ${TARGETS[0]} --output-csv ${TARGETS[1]} --dale True --rc True'
)

# Bound on rate of TM in the VB1-8 sequences with reverse complements included in donor set
env.Command(
    [os.path.join(OUTPUT_DIR, 'vb18_bound_rc.tex'), os.path.join(OUTPUT_DIR, 'vb18_bound_rc.csv')],
    [pmf_rate_vb18_vs_mouse_imgt_v, pmf_rate_gpt_vs_mouse_mock],
    'Rscript analysis/compute_bound.R --naive-rate ${SOURCES[0]} --gpt-fpr ${SOURCES[1]} --output-tex ${TARGETS[0]} --output-csv ${TARGETS[1]} --dale True --rc True'
)

## Probabilities of templated mutagenesis given different donor sets ##

# Probability of TM in the gpt sequences given TM from mock gpt set
prob_given_tm_from_gpt = env.Command(
    os.path.join(OUTPUT_DIR, 'prob_given_tm_from_gpt.csv'),
    [partis_gpt, os.path.join(DATA_DIR, 'reference_sets/gpt_mock_from_mouse.fasta')],
    'python analysis/compute_probs.py --input ${SOURCES[0]} --output $TARGET --references ${SOURCES[1]}'
)

# Probability of TM in the gpt sequences given TM from the 129S1 V gene set
prob_given_tm_from_v = env.Command(
    os.path.join(OUTPUT_DIR, 'prob_given_tm_from_v.csv'),
    [partis_gpt, os.path.join(DATA_DIR, 'reference_sets/mus_musculus_129S1_v_genes.fasta')],
    'python analysis/compute_probs.py --input ${SOURCES[0]} --output $TARGET --references ${SOURCES[1]}'
)

# Probability of TM in the gpt sequences given TM from mock gpt set with reverse complements included
prob_given_tm_from_gpt_rc = env.Command(
    os.path.join(OUTPUT_DIR, 'prob_given_tm_from_gpt_rc.csv'),
    [partis_gpt, os.path.join(DATA_DIR, 'reference_sets/gpt_mock_from_mouse.fasta')],
    'python analysis/compute_probs.py --input ${SOURCES[0]} --output $TARGET --references ${SOURCES[1]} --rc True'
)

# Probability of TM in the gpt sequences given TM from the 129S1 V gene set with reverse complements included
prob_given_tm_from_v_rc = env.Command(
    os.path.join(OUTPUT_DIR, 'prob_given_tm_from_v_rc.csv'),
    [partis_gpt, os.path.join(DATA_DIR, 'reference_sets/mus_musculus_129S1_v_genes.fasta')],
    'python analysis/compute_probs.py --input ${SOURCES[0]} --output $TARGET --references ${SOURCES[1]} --rc True'
)


# Per base probability of TM in the gpt sequences given TM from the gpt mock gene set based on mouse V genes
per_base_gpt_gpt_mock_from_mouse = env.Command(
    os.path.join(OUTPUT_DIR, 'per_base_gpt_gpt_mock_from_mouse.csv'),
    [partis_gpt, os.path.join(DATA_DIR, 'reference_sets/gpt_mock_from_mouse.fasta')],
    'python analysis/compute_prob_per_base.py --input ${SOURCES[0]} --output $TARGET --references ${SOURCES[1]}'
)

# Per base probability of TM in the gpt sequences given TM from the 129S1 V gene set
per_base_gpt_129S1 = env.Command(
    os.path.join(OUTPUT_DIR, 'per_base_gpt_129S1.csv'),
    [partis_gpt, os.path.join(DATA_DIR, 'reference_sets/mus_musculus_129S1_v_genes.fasta')],
    'python analysis/compute_prob_per_base.py --input ${SOURCES[0]} --output $TARGET --reference ${SOURCES[1]}'
)

## Plots ##

# Figure 1: Plot of fraction of mutations in the gpt gene explainable by TM from the mock gpt set or the mouse IMGT V set
env.Command(
    os.path.join(OUTPUT_DIR, 'gpt_v_vs_gpt.pdf'),
    [pmf_rate_gpt_vs_mouse_mock, pmf_rate_gpt_vs_mouse_imgt_v],
    'Rscript analysis/make_combined_rate_plot.R --input-1 ${SOURCES[0]} --input-2 ${SOURCES[1]} --output $TARGET --dale-method True'
)

# Figure 1 supplement version: Plot of fraction of mutations in the gpt gene explainable by TM from the mock gpt set or the mouse IMGT V set with reverse complements included in both
env.Command(
    os.path.join(OUTPUT_DIR, 'gpt_v_vs_gpt_rc.pdf'),
    [pmf_rate_gpt_vs_mouse_mock, pmf_rate_gpt_vs_mouse_imgt_v],
    'Rscript analysis/make_combined_rate_plot.R --input-1 ${SOURCES[0]} --input-2 ${SOURCES[1]} --output $TARGET --dale-method True --rc True'
)

# Figure 2: Average probability of the observed mutations under a TM model, either from the gpt set or from the muose 129S1 set
env.Command(
    os.path.join(OUTPUT_DIR, 'prob_given_gcv.pdf'),
    [prob_given_tm_from_gpt, prob_given_tm_from_v],
    'Rscript analysis/prob_given_gcv_plot.R --input-1 ${SOURCES[0]} --input-2 ${SOURCES[1]} --output $TARGET'
)

# Figure 2 supplemental version: Average probability of the observed mutations under a TM model, either from the gpt set or from the mouse 129S1 set with reverse complements included in both
env.Command(
    os.path.join(OUTPUT_DIR, 'prob_given_gcv_rc.pdf'),
    [prob_given_tm_from_gpt_rc, prob_given_tm_from_v_rc],
    'Rscript analysis/prob_given_gcv_plot.R --input-1 ${SOURCES[0]} --input-2 ${SOURCES[1]} --output $TARGET'
)

# Figure 3: Plot describing the per base observed vs. expected probabilities
env.Command(
    os.path.join(OUTPUT_DIR, 'per_base_obs_vs_exp.pdf'),
    [per_base_gpt_gpt_mock_from_mouse, per_base_gpt_129S1],
    'Rscript analysis/make_per_base_plot.R --input-1 ${SOURCES[0]} --input-2 ${SOURCES[1]} --output $TARGET'
)
