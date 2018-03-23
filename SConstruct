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

# compute the false positive rate for motif finder and poly motif finder on the gpt data
fpr_gpt = Command(
	['analysis/output/fpr.csv', 'analysis/output/fpr_poly.csv'],
	partis_gpt,
	'python analysis/compute_fpr_gpt.py -i $SOURCE -j ${TARGETS[0]} -k ${TARGETS[1]}')
# compute the false positive rate for motif finder and poly motif finder on the ebola data
fpr_ebola = Command(
	['analysis/output/fpr_ebola.csv', 'analysis/output/fpr_poly_ebola.csv'],
	partis_ebola,
	'python analysis/compute_fpr_ebola.py -i $SOURCE -j ${TARGETS[0]} -k ${TARGETS[1]}')

# make plots of the false positive rate for mf and pmf
Command(
	'analysis/output/fpr_gpt.pdf',
	fpr_gpt[0],
	'Rscript analysis/make_fpr_plot.R $SOURCE $TARGET')
Command(
	'analysis/output/fpr_poly_gpt.pdf',
	fpr_gpt[1],
	'Rscript analysis/make_fpr_plot.R $SOURCE $TARGET')
Command(
	'analysis/output/fpr_ebola.pdf',
	fpr_ebola[0],
	'Rscript analysis/make_fpr_plot.R $SOURCE $TARGET')
Command(
	'analysis/output/fpr_poly_ebola.pdf',
	fpr_ebola[1],
	'Rscript analysis/make_fpr_plot.R $SOURCE $TARGET')