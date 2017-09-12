import unittest
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from motif_finder import motif_finder
from motif_finder import longest_motif_finder

class testMF(unittest.TestCase):

    def setUp(self):
        pass

    def test_edge(self):
        # set up test data
        r1 = SeqRecord(Seq("GCC"), name="good_1")
        r2 = SeqRecord(Seq("CCC"), name="good_2")
        r3 = SeqRecord(Seq("CCG"), name="no_good_1")
        r4 = SeqRecord(Seq("ATG"), name="no_good_2")
        query = "ATGCCCCG"
        mut_map = {query: [4]}
        ref = [r1, r2, r3, r4]
        # run motif finder with window size 3
        mf_out = motif_finder(mut_map, ref, 3)
        # there should be two entries in the output, one for each mutation
        self.assertEqual(len(mf_out), 1)
        # the mutation at position 4 matches to good_1 and good_2
        self.assertEqual(mf_out[(query, 4)], set(["good_1", "good_2"]))

    def test_window_longer_than_query(self):
        # set up test data
        r1 = SeqRecord(Seq("CCCC"), name="r1")
        query = "ATGCCCC"
        mut_map = {query: [0, 4]}
        ref = [r1]
        # run motif finder with a window size of 10, longer than the
        # length of the query
        mf_out = motif_finder(mut_map, ref, 10)
        # we should have two entries, one for each mutation
        self.assertEqual(len(mf_out), 2)
        # both of the mutations should map to an empty set
        self.assertEqual(mf_out[(query, 0)], set())
        self.assertEqual(mf_out[(query, 4)], set())

    def test_multiple_queries(self):
        # set up test data
        r1 = SeqRecord(Seq("ATGCCCC"), name="r1")
        q1 = "CCCGGGGGG"
        q2 = "GGGGGGCCA"
        mut_map = {q1: [1, 4], q2: [5]}
        ref = [r1]
        # run motif finder
        mf_out = motif_finder(mut_map, ref, 3)
        self.assertEqual(len(mf_out), 3)
        self.assertEqual(mf_out[(q1, 1)], set(["r1"]))
        self.assertEqual(mf_out[(q1, 4)], set())
        self.assertEqual(mf_out[(q2, 5)], set(["r1"]))


class testLMF(unittest.TestCase):

    def setUp(self):
        pass

    def test_1(self):
        # set up the test data
        match = "ATGCCCC"
        query = "GGGGG" + match + "GGGGG"
        # r1 is the sequence that has a match to the query
        r1 = SeqRecord(Seq("AAAAAAAA" + match + "AAAAAAAA"),
                       name="r1")
        # r2 is a decoy sequence
        r2 = SeqRecord(Seq("AAAAAAAAAAAAAAAAAAAA"),
                       name="r2")
        # we want to find a match to the query around position 7
        mut_idx = 7
        mut_map = {query: [mut_idx]}
        ref = [r1, r2]
        seed_len = 3
        # run longest_motif_finder
        lmf_out = longest_motif_finder(mut_map, ref, seed_len)
        # the query centered at 7 matches to 1 position in the
        # reference
        self.assertEqual(len(lmf_out[(query, mut_idx)]), 1)
        # the match should start at position 8 and go for the
        # entire length of 'match'
        self.assertEqual(lmf_out[(query, mut_idx)],
                         set([("r1", str(r1.seq), 8, len(match))]))

    def test_overlapping(self):
        # set up the test data
        match = "A" * 10
        query = "G" * 10 + match + "G" * 10
        r1 = SeqRecord(Seq("C" * 10 + match + "C" * 10),
                       name="r1")
        mut_idx = 15
        mut_map = {query: [mut_idx]}
        ref = [r1]
        seed_len = 4
        # run longest_motif_finder
        lmf_out = longest_motif_finder(mut_map, ref, seed_len)
        # there should be a match to r1, starting at position 10, of
        # length len(match)
        self.assertEqual(("r1", str(r1.seq), 10, len(match)) in
                         lmf_out[(query, mut_idx)], True)


if __name__ == '__main__':
    unittest.main()
