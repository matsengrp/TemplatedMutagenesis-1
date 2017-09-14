import unittest
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from motif_finder import motif_finder
from motif_finder import longest_motif_finder
from motif_finder import seed_starts
from motif_finder import poly_motif_finder


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
        self.assertEqual(seed_starts(1, 1, 3, 10), (0, 1))
        # mutation at 5, window size 3, sequence of length 10 has
        # windows starting at 3, 4, 5 around the mutation
        self.assertEqual(seed_starts(5, 5, 3, 10), (3, 5))
        # mutation at 9, window size 3, sequence of length 10 has one
        # window starting at 7 around the mutation
        self.assertEqual(seed_starts(9, 9, 3, 10), (7, 7))

    def test_double_mutation(self):
        # mutations at 0 and 3, window size 4, sequence length 10 has
        # a window starting at 0 containing the two mutations
        self.assertEqual(seed_starts(0, 3, 4, 10), (0, 0))
        # mutations at 4 and 5, window size 4, sequence length 10 has
        # windows starting at 2, 3, and 4 containing the mutation
        self.assertEqual(seed_starts(4, 5, 4, 10), (2, 4))
        # mutations at 8 and 9, window size 4, sequence length 10 has
        # one window starting at 6 containing the mutations
        self.assertEqual(seed_starts(8, 9, 4, 10), (6, 6))


if __name__ == '__main__':
    unittest.main()
