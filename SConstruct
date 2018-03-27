#!/usr/bin/env scons

# use partis to call mutations on gpt sequences
partis_gpt = Command(
	   'run_partis/partis_output_gpt',
	   'data/yeap/presto_output',
	   'python run_partis/run_partis_gpt.py -i $SOURCE -o $TARGET')

# use partis to call mutations on ebola sequences
partis_ebola = Command(
	      'run_partis/partis_output_ebola',
	      'data/ebola/ebola_sequences_heavy.fasta',
	      'python run_partis/run_partis_ebola.py -i $SOURCE -o $TARGET')

# compute the false positive rate for motif finder and poly motif
# finder on the gpt sequences vs. gpt reference set
fpr_gpt_vs_gpt = Command(
	['analysis/output/fpr_gpt_gpt.csv',
         'analysis/output/fpr_poly_gpt_gpt.csv'],
	[partis_gpt,
         'data/reference_sets/gpt_132.fasta'],
	'python analysis/compute_fpr.py --input-directory ${SOURCES[0]} --output-csv ${TARGETS[0]} --output-csv-poly ${TARGETS[1]} --references ${SOURCES[1]}')

# compute the fpr for mf/pmf on the gpt sequences vs v reference set
fpr_gpt_vs_v = Command(
	['analysis/output/fpr_gpt_v.csv',
         'analysis/output/fpr_poly_gpt_v.csv'],
	[partis_gpt,
         'data/reference_sets/mus_musculus_129S1_v_genes.fasta'],
	'python analysis/compute_fpr.py --input-directory ${SOURCES[0]} --output-csv ${TARGETS[0]} --output-csv-poly ${TARGETS[1]} --references ${SOURCES[1]}')

# compute the fpr for mf/pmf on the ebola data vs. the v gene reference set
fpr_ebola = Command(
	['analysis/output/fpr_ebola.csv',
         'analysis/output/fpr_poly_ebola.csv'],
	[partis_ebola,
         'data/reference_sets/imgt_ighv_human.fasta'],
	'python analysis/compute_fpr.py --input-directory ${SOURCES[0]} --output-csv ${TARGETS[0]} --output-csv-poly ${TARGETS[1]} --references ${SOURCES[1]}')

# compute the average probability of mutation given gene conversion on the gpt data
prob_given_gcv = Command(
    'analysis/output/likelihood_gpt.csv',
    partis_gpt,
    'python analysis/compute_likelihoods.py -i $SOURCE -o $TARGET'
    )

# make plots of the false positive rate for mf and pmf
Command(
	'analysis/output/fpr_gpt.pdf',
	[fpr_gpt_vs_gpt[0], fpr_gpt_vs_v[0]],
	'Rscript analysis/make_fpr_plot.R --input-1 ${SOURCES[0]} --input-2 ${SOURCES[1]} --output $TARGET')
Command(
	'analysis/output/fpr_poly_gpt.pdf',
	[fpr_gpt_vs_gpt[1], fpr_gpt_vs_v[1]],
	'Rscript analysis/make_fpr_plot.R --input-1 ${SOURCES[0]} --input-2 ${SOURCES[1]} --output $TARGET')
