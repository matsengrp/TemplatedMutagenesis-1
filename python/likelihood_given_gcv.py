from process_partis import process_partis
from motif_finder import indexed_motif_finder, make_kmer_dictionary, extend_matches, n_alignments_per_mutation
from Bio import SeqIO
import numpy as np
import pandas as pd


def likelihood_given_gcv(partis_file, reference_file, k):
    """Finds the likelihood of mutations conditional on being due to gcv

    Keyword arguments:
    partis_file -- A partis csv describing the mutations.
    reference_file -- A fasta file with the reference sequences.
    k -- The minimum match length for gcv tracts.

    Returns: A data frame giving the probability of seeing the each
    observed mutation. Mutations are described by the name of the
    query sequence and the position of the mutation in that query
    sequence.

    """
    bases = ["A", "C", "G", "T"]
    mut_df = process_partis(partis_file)
    references = [r for r in SeqIO.parse(reference_file, "fasta")]
    kmer_dict = make_kmer_dictionary(references, k)
    # make a data frame containing all the mutations we didn"t see
    unobs_mut_rows = []
    for index, row in mut_df.iterrows():
        for b in bases:
            if b not in set([row["gl_base"], row["mutated_base"]]):
                r = row.copy()
                unseen_seq = list(r["mutated_seq"])
                unseen_seq[r["mutation_index"]] = b
                r["mutated_seq"] = "".join(unseen_seq)
                r["mutated_base"] = b
                unobs_mut_rows.append(r)
    unobs_mut_df = pd.DataFrame(unobs_mut_rows)

    # run motif finder on the observed and unobserved mutations
    motifs_obs = n_alignments_per_mutation(mut_df, kmer_dict, k)
    motifs_unobs = n_alignments_per_mutation(unobs_mut_df, kmer_dict, k)
    obs_and_unobs = pd.merge(motifs_obs, motifs_unobs, 
                             how="outer", 
                             on=["query_mutation_index", "query_name"],
                             validate="one_to_one")
    # get the probabilities of seeing the observed
    # mutations. n_alignments_x is the number of alignments for
    # observed mutations because motifs_obs was in the first position
    # in pd.merge
    def get_prob(row):
        n_obs = row["n_alignments_x"]
        n_unobs = row["n_alignments_y"]
        if n_obs + n_unobs == 0:
            return(np.nan)
        return n_obs / (n_obs + n_unobs + 0.)
    obs_and_unobs["prob"] = obs_and_unobs.apply(get_prob, axis=1)
    return(obs_and_unobs)
