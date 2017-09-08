import unittest
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from motif_finder import motif_finder


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


if __name__ == '__main__':
    unittest.main()
