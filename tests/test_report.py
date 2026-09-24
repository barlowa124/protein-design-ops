import unittest

from design_ops.report import identity, pairwise_identities, spearman


class ReportTests(unittest.TestCase):
    def test_identity(self):
        self.assertEqual(identity("ACGT", "ACGT"), 1.0)
        self.assertEqual(identity("ACGT", "AGGT"), 0.75)
        self.assertEqual(identity("AAAA", "TTTT"), 0.0)

    def test_pairwise(self):
        ids = pairwise_identities(["AAAA", "AATT", "AATT"])
        self.assertEqual(len(ids), 3)
        self.assertAlmostEqual(ids[0], 0.5)
        self.assertAlmostEqual(ids[2], 1.0)

    def test_spearman(self):
        self.assertAlmostEqual(spearman([1, 2, 3], [3, 2, 1]), -1.0)
        self.assertAlmostEqual(spearman([1, 2, 3], [1, 2, 3]), 1.0)


if __name__ == "__main__":
    unittest.main()
