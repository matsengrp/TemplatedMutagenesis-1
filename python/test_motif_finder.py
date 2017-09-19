import unittest
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from motif_finder import motif_finder
from motif_finder import longest_motif_finder
from motif_finder import seed_starts
from motif_finder import poly_motif_finder
from motif_finder import make_kmer_dictionary
from motif_finder import indexed_motif_finder
import pandas as pd


class testMotifFinder(unittest.TestCase):

    def setUp(self):
        pass

    def test_single(self):
        # set up test data
        match = "CCC"
        query = "GGG" + match + "GGG"
        r1 = SeqRecord(Seq("AAAAA" + match + "AAAAA"),
                       name="r1")
        mut_idx = 4
        mut_map = {query: [mut_idx]}
        ref = [r1]
        seed_len = 3
        # run motif_finder
        mf_out = motif_finder(mut_map, ref, seed_len)
        # there should be one match to r1
        self.assertEqual(mf_out[0].num_hits(), 1)
        self.assertEqual(mf_out[0].ref_idx, [6])
        self.assertEqual(len(mf_out), 1)

    def test_double(self):
        # set up test data
        match = "CCC"
        query = "GGG" + match + "GGG"
        r1 = SeqRecord(Seq("AAAAA" + match + "AAAAA" + match),
                       name="r1")
        mut_idx = 4
        mut_map = {query: [mut_idx]}
        ref = [r1]
        seed_len = 3
        # run motif_finder
        mf_out = motif_finder(mut_map, ref, seed_len)
        # there should be two matches to r1, at indices 6 and 14
        self.assertEqual(len(mf_out), 1)
        self.assertEqual(mf_out[0].num_hits(), 2)
        self.assertEqual(mf_out[0].ref_idx, [6, 14])


class testIndexedMotifFinder(unittest.TestCase):

    def setUp(self):
        pass

    def test_single(self):
        match = "CCC"
        query = "GGG" + match + "GGG"
        r1 = SeqRecord(Seq("AAAAA" + match + "AAAAA"),
                       name="r1")
        mut_idx = 4
        mut_map = {query: [mut_idx]}
        ref = [r1]
        k = 3
        kmer_dict = make_kmer_dictionary(ref, k)
        mf_out = indexed_motif_finder(mut_map, kmer_dict, k)
        self.assertEqual(mf_out.shape[0], 1)
        self.assertEqual(mf_out["reference_name"][0], "r1")
        self.assertEqual(mf_out["query_sequence"][0], query)
        self.assertEqual(mf_out["query_mutation_index"][0], mut_idx)

    def test_double(self):
        match = "CCC"
        query = "GGG" + match + "GGG"
        r1 = SeqRecord(Seq("AAAAA" + match + "AAAAA" + match),
                       name="r1")
        mut_idx = 4
        mut_map = {query: [mut_idx]}
        ref = [r1]
        k = 3
        kmer_dict = make_kmer_dictionary(ref, k)
        mf_out = indexed_motif_finder(mut_map, kmer_dict, k)
        self.assertEqual(mf_out.shape[0], 2)
        self.assertEqual(
            mf_out["query_sequence"].equals(pd.Series([query, query])), True)
        self.assertEqual(
            mf_out["query_mutation_index"].equals(pd.Series([mut_idx, mut_idx])), True)
        self.assertEqual(
            mf_out["reference_name"].equals(pd.Series(["r1", "r1"])), True)
        self.assertEqual(
            mf_out["reference_alignment"].equals(pd.Series([6, 14])), True)


class testPolyMotifFinder(unittest.TestCase):

    def setUp(self):
        pass

    def test_single(self):
        match = "CGC"
        query = "GGG" + match + "GGG"
        r1 = SeqRecord(Seq("AAAAA" + match + "AAAAA"),
                       name="r1")
        mut_idx = (3, 4)
        mut_map = {query: [mut_idx]}
        ref = [r1]
        seed_len = 3
        # run motif_finder
        pmf_out = poly_motif_finder(mut_map, ref, seed_len)
        # there should be one match to r1, at indices 5 and 6
        self.assertEqual(pmf_out[0].num_hits(), 1)
        self.assertEqual(pmf_out[0].ref_idx, [(5, 6)])
        self.assertEqual(len(pmf_out), 1)

    def test_double(self):
        match = "CCC"
        query = "GGG" + match + "GGG"
        r1 = SeqRecord(Seq("A" * 5 + match + "A" * 5 + match),
                       name="r1")
        mut_idx = (3, 5)
        mut_map = {query: [mut_idx]}
        ref = [r1]
        seed_len = 3
        pmf_out = poly_motif_finder(mut_map, ref, seed_len)
        self.assertEqual(len(pmf_out), 1)
        self.assertEqual(pmf_out[0].num_hits(), 2)
        self.assertEqual(pmf_out[0].ref_idx, [(5, 7), (13, 15)])


class testLongestMotifFinder(unittest.TestCase):

    def setUp(self):
        pass

    def test_single(self):
        # set up the test data
        match = "ATGCA"
        query = "G" * 5 + match + "G" * 5
        r1 = SeqRecord(Seq("A" * 10 + match + "A" * 10),
                       name="r1")
        mut_idx = 5
        mut_map = {query: [mut_idx]}
        ref = [r1]
        seed_len = 3
        # run longest_motif_finder
        lmf_out = longest_motif_finder(mut_map, ref, seed_len)
        self.assertEqual(len(lmf_out), 1)
        self.assertEqual(lmf_out[0].ref_name, ["r1"])
        self.assertEqual(lmf_out[0].ref_idx, [10])
        self.assertEqual(lmf_out[0].match_extent, [len(match)])

    def test_double(self):
        match = "ATGCA"
        query = "G" * 5 + match + "G" * 5
        r1 = SeqRecord(Seq("A" * 10 + match + "A" * 10 + match),
                       name="r1")
        mut_idx = 5
        mut_map = {query: [mut_idx]}
        ref = [r1]
        seed_len = 3
        # run longest_motif_finder
        lmf_out = longest_motif_finder(mut_map, ref, seed_len)
        self.assertEqual(len(lmf_out), 1)
        self.assertEqual(lmf_out[0].num_hits(), 2)
        self.assertEqual(lmf_out[0].ref_idx, [10, 25])
        self.assertEqual(lmf_out[0].match_extent, [len(match), len(match)])


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
