from process_partis import process_partis, unseen_mutations
from motif_finder import indexed_motif_finder, make_kmer_dictionary, extend_matches
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
    # read in the reference sequences and create the k-mer dictionary
    refs = [r for r in SeqIO.parse(reference_file, "fasta")]
    kmerdict = make_kmer_dictionary(refs, k)
    # get the observed mutations from the partis file
    mutations_obs = process_partis(partis_file)
    # the unobserved mutations
    mutations_unobs = unseen_mutations(partis_file)
    # run motif finder on the observed and unobserved mutations
    motifs_obs = indexed_motif_finder(mutations_obs, kmerdict, k)
    motifs_unobs = indexed_motif_finder(mutations_unobs, kmerdict, k)
    # extend all the matches and drop duplicates to get the unique
    # alignments for the observed and unobserved mutations
    extend_matches(motifs_obs)
    extend_matches(motifs_unobs)
    motifs_obs.drop_duplicates(inplace=True)
    motifs_unobs.drop_duplicates(inplace=True)
    # make a data frame containing all the (query, mutation_index)
    # pairs. For each pair we will find the number of matches for the
    # observed mutation and the unobserved mutations.
    query_and_idx = motifs_obs[["query_mutation_index", "query_name"]].drop_duplicates()
    rows = []
    for (q, i) in zip(query_and_idx["query_name"],
                      query_and_idx["query_mutation_index"]):
        # the number of matches is the number of times we have an
        # alignment in motifs_obs that is not the empty string (which
        # refers to no match)
        obs_templates = motifs_obs.loc[(motifs_obs.query_name == q) &
                                       (motifs_obs.query_mutation_index == i), :]
        n_obs = len([a for a in obs_templates["reference_name"] if a != ""])
        unobs_templates = motifs_unobs.loc[(motifs_unobs.query_name == q) &
                                     (motifs_unobs.query_mutation_index == i), :]
        n_unobs = len([a for a in unobs_templates["reference_name"] if a != ""])
        # if the top and bottom are both 0, the probability is not
        # defined (gene conversion was impossible because there were
        # no suitable templates)
        if (n_obs + n_unobs) == 0:
            prob = np.nan
        else:
            prob = float(n_obs) / float(n_obs + n_unobs)
        rows.append({"query_name": q,
                     "query_mutation_index": i,
                     "probability": prob})
    probdf = pd.DataFrame(rows)
    return(probdf)
