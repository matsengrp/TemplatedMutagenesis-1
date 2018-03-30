#!/usr/bin/env scons
# -*- python -*-

import os
DATA_DIR = 'data'
PARTIS = '/Users/juliefukuyama/GitHub/partis/bin/partis'


# use partis to call mutations on gpt sequences
partis_gpt = Command(
	   'run_partis/partis_output_gpt',
	   [os.path.join(DATA_DIR, 'yeap/presto_output'), PARTIS],
	   'python run_partis/run_partis_gpt.py --input-directory ${SOURCES[0]} --partis ${SOURCES[1]} --output-directory $TARGET')

# use partis to call mutations on ebola sequences
partis_ebola = Command(
	      'run_partis/partis_output_ebola',
	      [os.path.join(DATA_DIR, 'ebola/ebola_sequences_heavy.fasta'), PARTIS],
	      'python run_partis/run_partis_ebola.py --input-file ${SOURCES[0]} --partis ${SOURCES[1]} --output-directory $TARGET')

# compute the false positive rate for motif finder and poly motif
# finder on the gpt sequences vs. gpt reference set
fpr_gpt_vs_gpt = Command(
	['output/fpr_gpt_gpt.csv',
         'output/fpr_poly_gpt_gpt.csv'],
	[partis_gpt,
         os.path.join(DATA_DIR, 'reference_sets/gpt_132.fasta')],
	'python compute_fpr.py --input-directory ${SOURCES[0]} --output-csv ${TARGETS[0]} --output-csv-poly ${TARGETS[1]} --references ${SOURCES[1]}')

# compute the fpr for mf/pmf on the gpt sequences vs v reference set
fpr_gpt_vs_v = Command(
	['output/fpr_gpt_v.csv',
         'output/fpr_poly_gpt_v.csv'],
	[partis_gpt,
         os.path.join(DATA_DIR, 'reference_sets/mus_musculus_129S1_v_genes.fasta')],
	'python compute_fpr.py --input-directory ${SOURCES[0]} --output-csv ${TARGETS[0]} --output-csv-poly ${TARGETS[1]} --references ${SOURCES[1]}')

# compute the fpr for mf/pmf on the ebola data vs. the v gene reference set
fpr_ebola = Command(
	['output/fpr_ebola.csv',
         'output/fpr_poly_ebola.csv'],
	[partis_ebola,
         os.path.join(DATA_DIR, 'reference_sets/imgt_ighv_human.fasta')],
	'python compute_fpr.py --input-directory ${SOURCES[0]} --output-csv ${TARGETS[0]} --output-csv-poly ${TARGETS[1]} --references ${SOURCES[1]}')

# compute the probability of each mutation given gene conversion on the gpt data
prob_given_gcv_gpt = Command(
    'output/likelihood_gpt_gpt.csv',
    [partis_gpt,
     os.path.join(DATA_DIR, 'reference_sets/gpt_132.fasta')],
    'python compute_likelihoods.py --input ${SOURCES[0]} --output $TARGET --references ${SOURCES[1]}'
    )

prob_given_gcv_v = Command(
    'output/likelihood_gpt_v.csv',
    [partis_gpt,
     os.path.join(DATA_DIR, 'reference_sets/mus_musculus_129S1_v_genes.fasta')],
    'python compute_likelihoods.py --input ${SOURCES[0]} --output $TARGET --references ${SOURCES[1]}'
    )

# compute the probability of each mutation to each base on the gpt data
per_base_prob_gpt_gpt = Command(
    'output/per_base_gpt_gpt.csv',
    [partis_gpt,
     os.path.join(DATA_DIR, 'reference_sets/gpt_132.fasta')],
    'python compute_prob_per_base.py --input ${SOURCES[0]} --output $TARGET --references ${SOURCES[1]}'
)

per_base_prob_gpt_v = Command(
    'output/per_base_gpt_v.csv',
    [partis_gpt,
     os.path.join(DATA_DIR, 'reference_sets/mus_musculus_129S1_v_genes.fasta')],
    'python compute_prob_per_base.py --input ${SOURCES[0]} --output $TARGET --references ${SOURCES[1]}'
)

# make plots of the false positive rate for mf and pmf
Command(
	'output/fpr_gpt.pdf',
	[fpr_gpt_vs_gpt[0], fpr_gpt_vs_v[0]],
	'Rscript analysis/make_fpr_plot.R --input-1 ${SOURCES[0]} --input-2 ${SOURCES[1]} --output $TARGET')
Command(
	'output/fpr_poly_gpt.pdf',
	[fpr_gpt_vs_gpt[1], fpr_gpt_vs_v[1]],
	'Rscript analysis/make_fpr_plot.R --input-1 ${SOURCES[0]} --input-2 ${SOURCES[1]} --output $TARGET')

Command(
    'output/prob_given_gcv.pdf',
    [prob_given_gcv_gpt, prob_given_gcv_v],
    'Rscript analysis/prob_given_gcv_plot.R --input-1 ${SOURCES[0]} --input-2 ${SOURCES[1]} --output $TARGET')

Command(
    'output/per_base_obs_vs_exp.pdf',
    [per_base_prob_gpt_gpt, per_base_prob_gpt_v],
    'Rscript analysis/make_per_base_plot.R --input-1 ${SOURCES[0]} --input-2 ${SOURCES[1]} --output $TARGET')
