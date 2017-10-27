import unittest
from likelihood_given_gcv import likelihood_given_gcv
import numpy as np


class testLikelihood(unittest.TestCase):

    def setUp(self):
        pass

    def test_likelihood(self):
        partis_file = "/Users/juliefukuyama/GitHub/gcgcgc/test_data/likelihood_test_partis.csv"
        ref_file = "/Users/juliefukuyama/GitHub/gcgcgc/test_data/likelihood_test_reference.fasta"
        probs = likelihood_given_gcv(partis_file, ref_file, 3)
        prob_s1 = probs.loc[probs.query_name == "s1", "probability"]
        prob_s2 = probs.loc[probs.query_name == "s2", "probability"]
        prob_s3 = probs.loc[probs.query_name == "s3", "probability"]
        # s1 should have prob = 1/2, s2 prob = 0, s3 prob = nan
        self.assertEqual(prob_s1.item(), .5)
        self.assertEqual(prob_s2.item(), 0)
        self.assertEqual(np.isnan(prob_s3.item()), True)



if __name__ == '__main__':
    unittest.main()
