#!/usr/bin/env scons
# -*- python -*-

import os
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

# use partis to call mutations on VB1-8 sequences
partis_vb18 = env.Command(
    os.path.join(OUTPUT_DIR, 'partis/partis_output_vb18'),
    [os.path.join(DATA_DIR, 'yeap/presto_output'), PARTIS],
    'python analysis/run_partis_vb18.py --input-directory ${SOURCES[0]} --partis ${SOURCES[1]} --output-directory $TARGET')

# compute the pmf rate on the gpt sequences with the gpt reference set
pmf_rate_gpt_vs_gpt = env.Command(
    os.path.join(OUTPUT_DIR, 'pmf_rate_gpt_gpt.csv'),
    [partis_gpt,
     os.path.join(DATA_DIR, 'reference_sets/gpt_132.fasta')],
    'python analysis/compute_mf_rate.py --input-directory ${SOURCES[0]} --reference-fasta ${SOURCES[1]} --output-csv $TARGET')

# compute the pmf rate on the gpt sequences with the 129S1 v reference set
pmf_rate_gpt_vs_v = env.Command(
    os.path.join(OUTPUT_DIR, 'pmf_rate_gpt_v.csv'),
    [partis_gpt,
     os.path.join(DATA_DIR, 'reference_sets/mus_musculus_129S1_v_genes.fasta')],
    'python analysis/compute_mf_rate.py --input-directory ${SOURCES[0]} --reference-fasta ${SOURCES[1]} --output-csv $TARGET')

# compute the pmf rate for pmf on the gpt sequences with the imgt v reference set
pmf_rate_gpt_vs_imgt_v = env.Command(
    os.path.join(OUTPUT_DIR, 'pmf_rate_gpt_imgt_v.csv'),
    [partis_gpt,
     os.path.join(DATA_DIR, 'reference_sets/imgt_ighv_mouse.fasta')],
    'python analysis/compute_mf_rate.py --input-directory ${SOURCES[0]} --reference-fasta ${SOURCES[1]} --output-csv $TARGET')

# compute the pmf rate on the ebola data with the v gene reference set
pmf_rate_ebola = env.Command(
    os.path.join(OUTPUT_DIR, 'pmf_rate_ebola.csv'),
    [partis_ebola,
     os.path.join(DATA_DIR, 'reference_sets/imgt_ighv_human.fasta')],
    'python analysis/compute_mf_rate.py --input-directory ${SOURCES[0]} --reference-fasta ${SOURCES[1]} --output-csv $TARGET')

# compute the pmf rate for the vb18 sequences with the 129S1 v reference set
pmf_rate_vb18 = env.Command(
    os.path.join(OUTPUT_DIR, 'pmf_rate_vb18.csv'),
    [partis_vb18,
     os.path.join(DATA_DIR, 'reference_sets/mus_musculus_129S1_v_genes.fasta')],
    'python analysis/compute_mf_rate.py --input-directory ${SOURCES[0]} --reference-fasta ${SOURCES[1]} --output-csv $TARGET'
)

# compute the probability of each mutation given gene conversion on the gpt sequences with gpt templates
prob_given_gcv_gpt = env.Command(
    os.path.join(OUTPUT_DIR, 'prob_gpt_gpt.csv'),
    [partis_gpt,
     os.path.join(DATA_DIR, 'reference_sets/gpt_132.fasta')],
    'python analysis/compute_probs.py --input ${SOURCES[0]} --output $TARGET --references ${SOURCES[1]}'
    )

# compute the probability of each mutation given gene conversion on the gpt sequences with v gene templates
prob_given_gcv_v = env.Command(
    os.path.join(OUTPUT_DIR, 'prob_gpt_v.csv'),
    [partis_gpt,
     os.path.join(DATA_DIR, 'reference_sets/mus_musculus_129S1_v_genes.fasta')],
    'python analysis/compute_probs.py --input ${SOURCES[0]} --output $TARGET --references ${SOURCES[1]}'
    )

# compute the probability of each mutation given gene conversion on the gpt sequences with gpt templates and reverse complements included
prob_given_gcv_gpt_rc = env.Command(
    os.path.join(OUTPUT_DIR, 'prob_gpt_gpt_rc.csv'),
    [partis_gpt,
     os.path.join(DATA_DIR, 'reference_sets/gpt_132.fasta')],
    'python analysis/compute_probs.py --input ${SOURCES[0]} --output $TARGET --references ${SOURCES[1]} --rc True'
    )

# compute the probability of each mutation given gene conversion on the gpt sequences with v gene templates and reverse complements included
prob_given_gcv_v_rc = env.Command(
    os.path.join(OUTPUT_DIR, 'prob_gpt_v_rc.csv'),
    [partis_gpt,
     os.path.join(DATA_DIR, 'reference_sets/mus_musculus_129S1_v_genes.fasta')],
    'python analysis/compute_probs.py --input ${SOURCES[0]} --output $TARGET --references ${SOURCES[1]} --rc True'
    )

# compute the probability of each mutation to each base on the gpt sequences with gpt templates
per_base_prob_gpt_gpt = env.Command(
    os.path.join(OUTPUT_DIR, 'per_base_gpt_gpt.csv'),
    [partis_gpt,
     os.path.join(DATA_DIR, 'reference_sets/gpt_132.fasta')],
    'python analysis/compute_prob_per_base.py --input ${SOURCES[0]} --output $TARGET --references ${SOURCES[1]}'
)

# compute the probability of each mutation to each base on the gpt sequences with v templates
per_base_prob_gpt_v = env.Command(
    os.path.join(OUTPUT_DIR, 'per_base_gpt_v.csv'),
    [partis_gpt,
     os.path.join(DATA_DIR, 'reference_sets/mus_musculus_129S1_v_genes.fasta')],
    'python analysis/compute_prob_per_base.py --input ${SOURCES[0]} --output $TARGET --references ${SOURCES[1]}'
)

# compute the bound on the rate of templated mutagenesis on the ebola sequences
env.Command(
    [os.path.join(OUTPUT_DIR, 'ebola_bound_dale.tex'), os.path.join(OUTPUT_DIR, 'ebola_bound_dale.csv')],
    [pmf_rate_ebola, pmf_rate_gpt_vs_gpt],
    'Rscript analysis/compute_bound.R --naive-rate ${SOURCES[0]} --gpt-fpr ${SOURCES[1]} --output-tex ${TARGETS[0]} --output-csv ${TARGETS[1]} --dale True'
)
# compute the bound on the rate of templated mutagenesis on the ebola sequences with reverse complements included
env.Command(
    [os.path.join(OUTPUT_DIR, 'ebola_bound_dale_rc.tex'), os.path.join(OUTPUT_DIR, 'ebola_bound_dale_rc.csv')],
    [pmf_rate_ebola, pmf_rate_gpt_vs_gpt],
    'Rscript analysis/compute_bound.R --naive-rate ${SOURCES[0]} --gpt-fpr ${SOURCES[1]} --output-tex ${TARGETS[0]} --output-csv ${TARGETS[1]} --dale True --rc True'
)

# compute the bound on the rate of gene conversion on the vb18 sequences
env.Command(
    [os.path.join(OUTPUT_DIR, 'vb18_bound_dale.tex'), os.path.join(OUTPUT_DIR, 'vb18_bound_dale.csv')],
    [pmf_rate_vb18, pmf_rate_gpt_vs_gpt],
    'Rscript analysis/compute_bound.R --naive-rate ${SOURCES[0]} --gpt-fpr ${SOURCES[1]} --output-tex ${TARGETS[0]} --output-csv ${TARGETS[1]} --dale True'
)

# compute the bound on the rate of gene conversion on the vb18 sequences with reverse complements included
env.Command(
    [os.path.join(OUTPUT_DIR, 'vb18_bound_dale_rc.tex'), os.path.join(OUTPUT_DIR, 'vb18_bound_dale_rc.csv')],
    [pmf_rate_vb18, pmf_rate_gpt_vs_gpt],
    'Rscript analysis/compute_bound.R --naive-rate ${SOURCES[0]} --gpt-fpr ${SOURCES[1]} --output-tex ${TARGETS[0]} --output-csv ${TARGETS[1]} --dale True --rc True'
)

# make plot of pmf rate in gpt sequences with gpt templates vs gpt sequences with v gene templates
env.Command(
    os.path.join(OUTPUT_DIR, 'gpt_v_vs_gpt_dale.pdf'),
    [pmf_rate_gpt_vs_gpt, pmf_rate_gpt_vs_v],
    'Rscript analysis/make_combined_rate_plot.R --input-1 ${SOURCES[0]} --input-2 ${SOURCES[1]} --output $TARGETS --dale-method True'
)

# make plot of pmf rate in gpt sequences with gpt templates vs gpt sequences with v gene templates, reverse complements included
env.Command(
    os.path.join(OUTPUT_DIR, 'gpt_v_vs_gpt_dale_rc.pdf'),
    [pmf_rate_gpt_vs_gpt, pmf_rate_gpt_vs_v],
    'Rscript analysis/make_combined_rate_plot.R --input-1 ${SOURCES[0]} --input-2 ${SOURCES[1]} --output $TARGETS --dale-method True --rc True'
)

# plot for probability of mutation given gene conversion
env.Command(
    os.path.join(OUTPUT_DIR, 'prob_given_gcv.pdf'),
    [prob_given_gcv_gpt, prob_given_gcv_v],
    'Rscript analysis/prob_given_gcv_plot.R --input-1 ${SOURCES[0]} --input-2 ${SOURCES[1]} --output $TARGET'
)

# plot for the probability of mutation given gene conversion with reverse complements
env.Command(
    os.path.join(OUTPUT_DIR, 'prob_given_gcv_rc.pdf'),
    [prob_given_gcv_gpt_rc, prob_given_gcv_v_rc],
    'Rscript analysis/prob_given_gcv_plot.R --input-1 ${SOURCES[0]} --input-2 ${SOURCES[1]} --output $TARGET'
)

# plot describing the per base observed vs. expected probabilities
env.Command(
    os.path.join(OUTPUT_DIR, 'per_base_obs_vs_exp.pdf'),
    [per_base_prob_gpt_gpt, per_base_prob_gpt_v],
    'Rscript analysis/make_per_base_plot.R --input-1 ${SOURCES[0]} --input-2 ${SOURCES[1]} --output $TARGET'
)

# stouffer plots
env.Command(
    [os.path.join(OUTPUT_DIR, 'two_betas.pdf'), os.path.join(OUTPUT_DIR, 'stouffer_distributions.pdf')],
    [],
    'Rscript analysis/stouffer_simulations.R --beta-output ${TARGETS[0]} --distribution-output ${TARGETS[1]}'
)

