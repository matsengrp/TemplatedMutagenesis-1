import unittest
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from motif_finder import seed_starts, make_kmer_dictionary, indexed_motif_finder, extend_matches, hit_fraction, n_alignments_per_mutation
from process_partis import process_partis
import pandas as pd


class testMotifFinder(unittest.TestCase):

    def setUp(self):
        pass

    def test_imf(self):
        partis_file = "/Users/juliefukuyama/GitHub/gcgcgc/test_data/partis_test.csv"
        mut_df = process_partis(partis_file)
        r1 = SeqRecord("ATA", name="r1")
        r2 = SeqRecord("CAA", name="r2")
        kmer_dict = make_kmer_dictionary([r1, r2], k=2)
        mut_df = process_partis(partis_file)
        imf = indexed_motif_finder(mut_df, kmer_dict, k=2)
        # the partis file has a naive sequence AAAAAAAG and a mutated
        # sequence AAAAAAAA, so we should have one hit to r2
        self.assertEqual(imf.shape[0], 1)
        self.assertEqual(imf["reference_name"][0], "r2")
        self.assertEqual(imf["reference_alignment"][0], 2)
        self.assertEqual(imf["query_name"][0], "s1")
        self.assertEqual(imf["query_mutation_index"][0], 7)

    def test_alignments_per_mutation(self):
        partis_file = "/Users/juliefukuyama/GitHub/gcgcgc/test_data/partis_test.csv"
        mut_df = process_partis(partis_file)
        r1 = SeqRecord("ATA", name="r1")
        r2 = SeqRecord("AAA", name="r2")
        kmer_dict = make_kmer_dictionary([r1, r2], k=2)
        mut_df = process_partis(partis_file)
        nalign = n_alignments_per_mutation(mut_df, kmer_dict, k=2)
        # there is one mutation in the partis file, and it can be
        # explained by two alignments from r2
        self.assertEqual(nalign.shape[0], 1)
        self.assertEqual(nalign["n_alignments"][0], 2)
        self.assertEqual(nalign["query_mutation_index"][0], 7)
        self.assertEqual(nalign["query_name"][0], "s1")

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
        r1 = SeqRecord(Seq("ATA"), name="r1")
        d1 = make_kmer_dictionary([r1], 2)
        d2 = make_kmer_dictionary([r1], 3)
        d3 = make_kmer_dictionary([r1], 4)
        self.assertEqual(len(d1.keys()), 2)
        self.assertEqual("AT" in d1.keys(), True)
        self.assertEqual("TA" in d1.keys(), True)
        self.assertEqual(d1["AT"], set([(r1, 0)]))
        self.assertEqual(d1["TA"], set([(r1, 1)]))
        self.assertEqual(len(d2.keys()), 1)
        self.assertEqual("ATA" in d2.keys(), True)
        self.assertEqual(d2["ATA"], set([(r1, 0)]))
        self.assertEqual(len(d3.keys()), 0)

    def test_multiple_refs(self):
        r1 = SeqRecord("ATA", name="r1")
        r2 = SeqRecord("CTA", name="r2")
        d1 = make_kmer_dictionary([r1, r2], 3)
        d2 = make_kmer_dictionary([r1, r2], 2)
        self.assertEqual(len(d1.keys()), 2)
        self.assertEqual(d1["ATA"], set([(r1, 0)]))
        self.assertEqual(d1["CTA"], set([(r2, 0)]))
        self.assertEqual(len(d2.keys()), 3)
        self.assertEqual(d2["AT"], set([(r1, 0)]))
        self.assertEqual(d2["TA"], set([(r1, 1), (r2, 1)]))
        self.assertEqual(d2["CT"], set([(r2, 0)]))


if __name__ == '__main__':
    unittest.main()
