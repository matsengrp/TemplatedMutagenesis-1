import unittest
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
from motif_finder import seed_starts, make_kmer_dictionary, indexed_motif_finder, extend_matches, hit_fraction, n_alignments_per_mutation, likelihood_given_gcv
import pandas as pd
import numpy as np


class testMotifFinder(unittest.TestCase):

    def setUp(self):
        pass

    def test_imf(self):
        mut_df = pd.DataFrame([{
            "mutated_seq": "AAAAAAAA",
            "naive_seq": "AAAAAAAG",
            "mutated_seq_id": "s1",
            "mutation_index": 7,
            "gl_base": "G",
            "mutated_base": "A"
        }])
        r1 = SeqRecord("ATA", name="r1")
        r2 = SeqRecord("CAA", name="r2")
        kmer_dict = make_kmer_dictionary([r1, r2], k=2)
        imf = indexed_motif_finder(mut_df, kmer_dict, k=2)
        # the partis file has a naive sequence AAAAAAAG and a mutated
        # sequence AAAAAAAA, so we should have one hit to r2
        self.assertEqual(imf.shape[0], 1)
        self.assertEqual(imf["reference_name"][0], "r2")
        self.assertEqual(imf["reference_alignment"][0], 2)
        self.assertEqual(imf["query_name"][0], "s1")
        self.assertEqual(imf["query_mutation_index"][0], 7)

    def test_no_alignments(self):
        mut_df = pd.DataFrame([{
            "mutated_seq": "AAAAAAAA",
            "naive_seq": "AAAAAAAG",
            "mutated_seq_id": "s1",
            "mutation_index": 7,
            "gl_base": "G",
            "mutated_base": "A"
        }])
        r1 = SeqRecord("TTTT", name="r1")
        r2 = SeqRecord("TCTC", name="r2")
        kmer_dict = make_kmer_dictionary([r1, r2], k=3)
        nalign = n_alignments_per_mutation(mut_df, kmer_dict, k=3)
        # there are no AAA sequences in the references (r1 and r2), so
        # we should get no alignments
        self.assertEqual(nalign.loc[nalign.query_name == "s1", "n_alignments"].item(), 0)

    def test_one_alignment(self):
        mut_df = pd.DataFrame([{
            "mutated_seq": "AAAAAAAA",
            "naive_seq": "AAAAAAAG",
            "mutated_seq_id": "s1",
            "mutation_index": 7,
            "gl_base": "G",
            "mutated_base": "A"
        }])
        r1 = SeqRecord("TTTT", name="r1")
        r2 = SeqRecord("TAAA", name="r2")
        kmer_dict = make_kmer_dictionary([r1, r2], k=3)
        nalign = n_alignments_per_mutation(mut_df, kmer_dict, k=3)
        # there is one AAA sequence in r1 that could serve as a
        # template for the mutation
        self.assertEqual(nalign.loc[nalign.query_name == "s1", "n_alignments"].item(), 1)

    def test_two_alignments_different_references(self):
        mut_df = pd.DataFrame([{
            "mutated_seq": "AAAAAAAA",
            "naive_seq": "AAAAAAAG",
            "mutated_seq_id": "s1",
            "mutation_index": 7,
            "gl_base": "G",
            "mutated_base": "A"
        }])
        r1 = SeqRecord("AAAT", name="r1")
        r2 = SeqRecord("TAAA", name="r2")
        kmer_dict = make_kmer_dictionary([r1, r2], k=3)
        nalign = n_alignments_per_mutation(mut_df, kmer_dict, k=3)
        # in this setup, there is one template for the mutation in r2
        self.assertEqual(nalign.loc[nalign.query_name == "s1", "n_alignments"].item(), 2)

    def test_two_alignments_same_reference(self):
        mut_df = pd.DataFrame([{
            "mutated_seq": "AAAAAAAA",
            "naive_seq": "AAAAAAAG",
            "mutated_seq_id": "s1",
            "mutation_index": 7,
            "gl_base": "G",
            "mutated_base": "A"
        }])
        r1 = SeqRecord("TTTT", name="r1")
        r2 = SeqRecord("TAAAGAAA", name="r2")
        kmer_dict = make_kmer_dictionary([r1, r2], k=3)
        nalign = n_alignments_per_mutation(mut_df, kmer_dict, k=3)
        # two templates for the mutation, both from r2
        self.assertEqual(nalign.loc[nalign.query_name == "s1", "n_alignments"].item(), 2)

    def test_two_alignments_overlapping(self):
        mut_df = pd.DataFrame([{
            "mutated_seq": "AAAAAAAA",
            "naive_seq": "AAAAAAAG",
            "mutated_seq_id": "s1",
            "mutation_index": 7,
            "gl_base": "G",
            "mutated_base": "A"
        }])
        r1 = SeqRecord("TTTT", name="r1")
        r2 = SeqRecord("TAAAA", name="r2")
        kmer_dict = make_kmer_dictionary([r1, r2], k=3)
        nalign = n_alignments_per_mutation(mut_df, kmer_dict, k=3)
        # two overlapping templates for the mutation, both in r2
        self.assertEqual(nalign.loc[nalign.query_name == "s1", "n_alignments"].item(), 2)

    def test_imf_poly(self):
        mut_df = pd.DataFrame([{
            "mutated_seq": "AAAAAAAA",
            "naive_seq": "AAAAAGAG",
            "mutated_seq_id": "s1",
            "mutation_index": (5, 7),
            "gl_base": "GG",
            "mutated_base": "AA"
        }])
        r1 = SeqRecord("AAAA", name="r1")
        r2 = SeqRecord("AGAG", name="r2")
        kmer_dict = make_kmer_dictionary([r1, r2], k=4)
        imf = indexed_motif_finder(mut_df, kmer_dict, k=4)
        # the partis file has a naive sequence AAAAAGAG and a mutated
        # sequence AAAAAAAA, so we should have one hit to r1
        self.assertEqual(imf.shape[0], 1)
        self.assertEqual(imf["reference_name"][0], "r1")
        self.assertEqual(imf["reference_alignment"][0], 1)
        self.assertEqual(imf["query_name"][0], "s1")
        self.assertEqual(imf["query_mutation_index"][0], (5,7))


    def test_likelihood(self):
        # set up
        partis_file = "../test_data/likelihood_test_partis.csv"
        ref_file = "../test_data/likelihood_test_reference.fasta"
        references = [r for r in SeqIO.parse(ref_file, "fasta")]
        refdict = make_kmer_dictionary(references, 3)
        probs = likelihood_given_gcv(partis_file, refdict, 3)
        prob_s1 = probs.loc[probs.query_name == "s1", "prob"]
        prob_s2 = probs.loc[probs.query_name == "s2", "prob"]
        prob_s3 = probs.loc[probs.query_name == "s3", "prob"]
        # s1 should have prob = 1/2, s2 prob = 0, s3 prob = nan
        self.assertEqual(prob_s1.item(), .5)
        self.assertEqual(prob_s2.item(), 0)
        self.assertEqual(np.isnan(prob_s3.item()), True)

    def test_hit_fraction(self):
        mut_df = pd.DataFrame([{
            "mutated_seq": "AAAAAAAA",
            "naive_seq": "AAAAAAAG",
            "mutated_seq_id": "s1",
            "mutation_index": 7,
            "gl_base": "G",
            "mutated_base": "A"
        }, {
            "mutated_seq": "TCTAAAAA",
            "naive_seq": "ACTAAAAA",
            "mutated_seq_id": "s2",
            "mutation_index": 0,
            "gl_base": "A",
            "mutated_base": "T"
        }])
        r1 = SeqRecord("TCTC", name="r1")
        kmer_dict = make_kmer_dictionary([r1], k=3)
        imf = indexed_motif_finder(mut_df, kmer_dict, k=3)
        # the first mutation has no templates in the reference and the
        # second one does, so we should get a hit fraction of .5.
        self.assertEqual(hit_fraction(imf), .5)

    def test_mutation_overlap(self):
        # here we want to test that indexed_motif_finder deals with
        # multiple mutations in a window correctly: if there are two
        # mutations next to each other and we want to know if one of
        # them is templated, we should be searching for templates that
        # contain only one of the mutations, not both
        mut_df = pd.DataFrame([{
            "mutated_seq": "GTGGG",
            "naive_seq": "AAGGG",
            "mutated_seq_id": "s1",
            "mutation_index": 0,
            "gl_base": "A",
            "mutated_base": "G"
        }])
        r1 = SeqRecord("GAG", name="r1")
        kmer_dict = make_kmer_dictionary([r1], k=3)
        imf = indexed_motif_finder(mut_df, kmer_dict, k=3)
        # there's only one mutation, and it should have a template, so
        # hit_fraction = 1
        self.assertEqual(hit_fraction(imf), 1)


class testSeedStarts(unittest.TestCase):

    def setUp(self):
        pass

    def test_single_mutation(self):
        # mutation at 1, window size 3, sequence length 10 has windows
        # starting at 0 and 1 around the mutation
        self.assertEqual(seed_starts(1, 3, 10), (0, 1))
        # mutation at 5, window size 3, sequence of length 10 has
        # windows starting at 3, 4, 5 around the mutation
        self.assertEqual(seed_starts(5, 3, 10), (3, 5))
        # mutation at 9, window size 3, sequence of length 10 has one
        # window starting at 7 around the mutation
        self.assertEqual(seed_starts(9, 3, 10), (7, 7))

    def test_double_mutation(self):
        # mutations at 0 and 3, window size 4, sequence length 10 has
        # a window starting at 0 containing the two mutations
        self.assertEqual(seed_starts((0, 3), 4, 10), (0, 0))
        # mutations at 4 and 5, window size 4, sequence length 10 has
        # windows starting at 2, 3, and 4 containing the mutation
        self.assertEqual(seed_starts((4, 5), 4, 10), (2, 4))
        # mutations at 8 and 9, window size 4, sequence length 10 has
        # one window starting at 6 containing the mutations
        self.assertEqual(seed_starts((8, 9), 4, 10), (6, 6))


class testKmerDict(unittest.TestCase):
    def setUp(self):
        pass

    def test_single_ref(self):
        # make k-mer dictionaries with k = 2,3,4 for the sequence ATA
        r1 = SeqRecord(Seq("ATA"), name="r1")
        d2 = make_kmer_dictionary([r1], 2)
        d3 = make_kmer_dictionary([r1], 3)
        d4 = make_kmer_dictionary([r1], 4)
        # the 2-mer dictionary should have AT and TA
        self.assertEqual(len(d2.keys()), 2)
        self.assertEqual("AT" in d2.keys(), True)
        self.assertEqual("TA" in d2.keys(), True)
        # the AT 2-mer was in the r1 sequence at position 0
        self.assertEqual(d2["AT"], set([(r1, 0)]))
        # the TA 2-mer was in the r1 sequence at position 1
        self.assertEqual(d2["TA"], set([(r1, 1)]))
        # the 3-mer dictionary should contain just ATA
        self.assertEqual(len(d3.keys()), 1)
        self.assertEqual("ATA" in d3.keys(), True)
        # the ATA 3-mer occurred in r1 at position 0
        self.assertEqual(d3["ATA"], set([(r1, 0)]))
        # there are no 4-mers in a sequence of size 3
        self.assertEqual(d4, {})

    def test_multiple_refs(self):
        # test making a dictionary out of multiple reference sequences
        r1 = SeqRecord("ATA", name="r1")
        r2 = SeqRecord("CTA", name="r2")
        # dictionary of 3-mers
        d3 = make_kmer_dictionary([r1, r2], 3)
        # dictionary of 2-mers
        d2 = make_kmer_dictionary([r1, r2], 2)
        # the 3-mer dictionary should have ATA at position 0 in r1 and
        # CTA in position 0 in r2
        self.assertEqual(len(d3.keys()), 2)
        self.assertEqual(d3["ATA"], set([(r1, 0)]))
        self.assertEqual(d3["CTA"], set([(r2, 0)]))
        # the 2-mer dictionary should have AT at position 0 in r1, TA
        # at position 1 in both r1 and r2, and CT at position 0 in r2
        self.assertEqual(len(d2.keys()), 3)
        self.assertEqual(d2["AT"], set([(r1, 0)]))
        self.assertEqual(d2["TA"], set([(r1, 1), (r2, 1)]))
        self.assertEqual(d2["CT"], set([(r2, 0)]))


if __name__ == '__main__':
    unittest.main()
