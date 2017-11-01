import unittest
from process_partis import get_pairs, unseen_mutations, process_partis
from likelihood_given_gcv import likelihood_given_gcv
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

class testPP(unittest.TestCase):

    def setUp(self):
        pass

    def test_process_partis(self):
        # run partis on the test data
        partis_file = "/Users/juliefukuyama/GitHub/gcgcgc/test_data/partis_test.csv"
        mut_df = process_partis(partis_file)
        # there is one mutation in the test at position 7 with gl base
        # G and mutated base A
        self.assertEqual(mut_df.shape[0], 1)
        self.assertEqual(mut_df["mutated_seq"][0], "AAAAAAAA")
        self.assertEqual(mut_df["mutated_seq_id"][0], "s1")
        self.assertEqual(mut_df["mutation_index"][0], 7)
        self.assertEqual(mut_df["gl_base"][0], "G")
        self.assertEqual(mut_df["mutated_base"][0], "A")

    def test_get_pairs(self):
        # set up the test data
        mut_indices = [0, 1, 8, 9, 10, 15, 18, 21, 23]
        # run get_pairs
        pairs = get_pairs(mut_indices, min_spacing=2)
        # check that 
        self.assertEqual((0, 1) in pairs, True)
        self.assertEqual((8, 9) in pairs, True)
        self.assertEqual((8, 10) in pairs, True)
        self.assertEqual((9, 10) in pairs, True)
        self.assertEqual((21, 23) in pairs, True)
        self.assertEqual(len(pairs), 5)


class testUnseenMutations(unittest.TestCase):

    def setUp(self):
        pass

    def test_correct_mutations(self):
        partis_file = "/Users/juliefukuyama/GitHub/gcgcgc/test_data/partis_test.csv"
        mm = unseen_mutations(partis_file)
        unseen_bases = set([k[mm[k][0]] for k in mm.keys()])
        unseen_seq_1 = Seq("AAAAAAAC")
        unseen_seq_2 = Seq("AAAAAAAT")
        # there is one mutated base, so there are two corresponding
        # unseen mutations
        self.assertEqual(len(mm), 2)
        # the unseen bases should be T and C
        self.assertEqual(unseen_bases, set(["T", "C"]))
        # the sequences should be AAAAAAAC and AAAAAAAT
        self.assertEqual(set([k.seq for k in mm.keys()]),
                         set([unseen_seq_1, unseen_seq_2]))
        # the ids for both of the sequences should be s1
        self.assertEqual(set([k.id for k in mm.keys()]),
                         set(["s1"]))

    def test_correct_mutations_indels(self):
        partis_file = "/Users/juliefukuyama/GitHub/gcgcgc/test_data/partis_test_indels.csv"
        mm = unseen_mutations(partis_file)
        unseen_bases = set([k[mm[k][0]] for k in mm.keys()])
        unseen_seq_1 = Seq("AAAAAAAC")
        unseen_seq_2 = Seq("AAAAAAAT")
        # there is one mutated base, so there are two corresponding
        # unseen mutations
        self.assertEqual(len(mm), 2)
        # the unseen bases should be T and C
        self.assertEqual(unseen_bases, set(["T", "C"]))
        # the sequences should be AAAAAAAC and AAAAAAAT
        self.assertEqual(set([k.seq for k in mm.keys()]),
                         set([unseen_seq_1, unseen_seq_2]))
        # the ids for both of the sequences should be s1
        self.assertEqual(set([k.id for k in mm.keys()]),
                         set(["s1"]))


if __name__ == '__main__':
    unittest.main()
