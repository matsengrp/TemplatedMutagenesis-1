#!/usr/bin/env scons
# -*- python -*-

import os
import sys
print(sys.path)
env = Environment(ENV = os.environ)
DATA_DIR = 'data'
OUTPUT_DIR = 'output'
PARTIS = 'partis/bin/partis'

if not os.path.exists(OUTPUT_DIR):
    os.mkdir(OUTPUT_DIR)

# use partis to call mutations on gpt sequences
partis_gpt = env.Command(
    os.path.join(OUTPUT_DIR, 'partis/partis_output_gpt'),
    [os.path.join(DATA_DIR, 'yeap/presto_output'), PARTIS],
    'python analysis/run_partis_gpt.py --input-directory ${SOURCES[0]} --partis ${SOURCES[1]} --output-directory $TARGET')

# use partis to call mutations on ebola sequences
partis_ebola = env.Command(
    os.path.join(OUTPUT_DIR, 'partis/partis_output_ebola'),
    [os.path.join(DATA_DIR, 'ebola/ebola_sequences_heavy.fasta'), PARTIS],
    'python analysis/run_partis_ebola.py --input-file ${SOURCES[0]} --partis ${SOURCES[1]} --output-directory $TARGET')

# compute the false positive rate for motif finder and poly motif
# finder on the gpt sequences vs. gpt reference set
fpr_gpt_vs_gpt = env.Command(
    os.path.join(OUTPUT_DIR, 'fpr_gpt_gpt.csv'),
    [partis_gpt,
     os.path.join(DATA_DIR, 'reference_sets/gpt_132.fasta')],
    'python analysis/compute_mf_rate.py --input-directory ${SOURCES[0]} --reference-fasta ${SOURCES[1]} --output-csv $TARGET')

# compute the fpr for mf/pmf on the gpt sequences vs 129S1 v reference set
fpr_gpt_vs_v = env.Command(
    os.path.join(OUTPUT_DIR, 'fpr_gpt_v.csv'),
    [partis_gpt,
     os.path.join(DATA_DIR, 'reference_sets/mus_musculus_129S1_v_genes.fasta')],
    'python analysis/compute_mf_rate.py --input-directory ${SOURCES[0]} --reference-fasta ${SOURCES[1]} --output-csv $TARGET')

# compute the fpr for mf/pmf on the gpt sequences vs imgt v reference set
fpr_gpt_vs_imgt_v = env.Command(
    os.path.join(OUTPUT_DIR, 'fpr_gpt_imgt_v.csv'),
    [partis_gpt,
     os.path.join(DATA_DIR, 'reference_sets/imgt_ighv_mouse.fasta')],
    'python analysis/compute_mf_rate.py --input-directory ${SOURCES[0]} --reference-fasta ${SOURCES[1]} --output-csv $TARGET')

# compute the fpr for mf/pmf on the ebola data vs. the v gene reference set
fpr_ebola = env.Command(
    os.path.join(OUTPUT_DIR, 'fpr_ebola.csv'),
    [partis_ebola,
     os.path.join(DATA_DIR, 'reference_sets/imgt_ighv_human.fasta')],
    'python analysis/compute_mf_rate.py --input-directory ${SOURCES[0]} --reference-fasta ${SOURCES[1]} --output-csv $TARGET')

# compute the probability of each mutation given gene conversion on the gpt data
prob_given_gcv_gpt = env.Command(
    os.path.join(OUTPUT_DIR, 'prob_gpt_gpt.csv'),
    [partis_gpt,
     os.path.join(DATA_DIR, 'reference_sets/gpt_132.fasta')],
    'python analysis/compute_probs.py --input ${SOURCES[0]} --output $TARGET --references ${SOURCES[1]}'
    )

prob_given_gcv_v = env.Command(
    os.path.join(OUTPUT_DIR, 'prob_gpt_v.csv'),
    [partis_gpt,
     os.path.join(DATA_DIR, 'reference_sets/mus_musculus_129S1_v_genes.fasta')],
    'python analysis/compute_probs.py --input ${SOURCES[0]} --output $TARGET --references ${SOURCES[1]}'
    )

# compute the probability of each mutation to each base on the gpt data
per_base_prob_gpt_gpt = env.Command(
    os.path.join(OUTPUT_DIR, 'per_base_gpt_gpt.csv'),
    [partis_gpt,
     os.path.join(DATA_DIR, 'reference_sets/gpt_132.fasta')],
    'python analysis/compute_prob_per_base.py --input ${SOURCES[0]} --output $TARGET --references ${SOURCES[1]}'
)

per_base_prob_gpt_v = env.Command(
    os.path.join(OUTPUT_DIR, 'per_base_gpt_v.csv'),
    [partis_gpt,
     os.path.join(DATA_DIR, 'reference_sets/mus_musculus_129S1_v_genes.fasta')],
    'python analysis/compute_prob_per_base.py --input ${SOURCES[0]} --output $TARGET --references ${SOURCES[1]}'
)

# compute the bound on the rate of gene conversion on the ebola data
ebola_bound = env.Command(
    [os.path.join(OUTPUT_DIR, 'ebola_bound.tex'), os.path.join(OUTPUT_DIR, 'ebola_bound.csv')],
    [fpr_ebola, fpr_gpt_vs_gpt],
    'Rscript analysis/compute_bound.R --ebola-rate ${SOURCES[0]} --gpt-fpr ${SOURCES[1]} --output-tex ${TARGETS[0]} --output-csv ${TARGETS[1]}'
)

# make plots of the false positive rate for mf and pmf
fpr_gpt_plot = env.Command(
    [os.path.join(OUTPUT_DIR, 'fpr_gpt.pdf'), os.path.join(OUTPUT_DIR, 'fpr_poly_gpt.pdf')],
    fpr_gpt_vs_gpt,
    'Rscript analysis/make_rate_plot.R --input $SOURCES --output-mf ${TARGETS[0]} --output-pmf ${TARGETS[1]}'
)

fpr_gpt_v_plot = env.Command(
    [os.path.join(OUTPUT_DIR, 'fpr_gpt_v.pdf'), os.path.join(OUTPUT_DIR, 'fpr_poly_gpt_v.pdf')],
    fpr_gpt_vs_v[0],
    'Rscript analysis/make_rate_plot.R --input $SOURCES --output-mf ${TARGETS[0]} --output-pmf ${TARGETS[1]}'
)

fpr_gpt_imgt_v_plot = env.Command(
    [os.path.join(OUTPUT_DIR, 'fpr_gpt_imgt_v.pdf'), os.path.join(OUTPUT_DIR, 'fpr_poly_gpt_imgt_v.pdf')],
    fpr_gpt_vs_imgt_v[0],
    'Rscript analysis/make_rate_plot.R --input $SOURCES --output-mf ${TARGETS[0]} --output-pmf ${TARGETS[1]}'
)

prob_given_gcv_plot = env.Command(
    os.path.join(OUTPUT_DIR, 'prob_given_gcv.pdf'),
    [prob_given_gcv_gpt, prob_given_gcv_v],
    'Rscript analysis/prob_given_gcv_plot.R --input-1 ${SOURCES[0]} --input-2 ${SOURCES[1]} --output $TARGET'
)

per_base_plot = env.Command(
    os.path.join(OUTPUT_DIR, 'per_base_obs_vs_exp.pdf'),
    [per_base_prob_gpt_gpt, per_base_prob_gpt_v],
    'Rscript analysis/make_per_base_plot.R --input-1 ${SOURCES[0]} --input-2 ${SOURCES[1]} --output $TARGET'
)

# make phylogenetic tree plots
tree_plots = env.Command(
    os.path.join(OUTPUT_DIR, 'gene_tree_plots.pdf'),
    [],
    'Rscript analysis/tree_plots.R --tree-output $TARGET'
)
