import unittest
from process_partis import get_pairs


class testPP(unittest.TestCase):

    def setUp(self):
        pass

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


if __name__ == '__main__':
    unittest.main()
