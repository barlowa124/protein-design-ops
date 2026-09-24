import tempfile
import unittest

from design_ops.backbone import parse_chain


def _ca(serial, resseq, resname, chain="A"):
    return (
        f"ATOM  {serial:5d}  CA  {resname:>3s} {chain}{resseq:4d}    "
        f"{0.0:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00 50.00           C\n"
    )


def _pdb(lines):
    tmp = tempfile.NamedTemporaryFile("w", suffix=".pdb", delete=False)
    tmp.writelines(lines)
    tmp.close()
    return tmp.name


class BackboneTests(unittest.TestCase):
    def test_parse_chain_sequence(self):
        pdb = _pdb([
            _ca(1, 1, "ASN"), _ca(2, 2, "LEU"), _ca(3, 3, "TYR"),
            _ca(4, 1, "GLY", "B"),  # other chain ignored
        ])
        info = parse_chain(pdb, "A")
        self.assertEqual(info["native_seq"], "NLY")
        self.assertEqual(info["resseqs"], [1, 2, 3])

    def test_missing_chain_raises(self):
        pdb = _pdb([_ca(1, 1, "ALA")])
        with self.assertRaises(ValueError):
            parse_chain(pdb, "Z")

    def test_real_fixture(self):
        info = parse_chain(
            "../dockops/tests/fixtures/1L2Y.pdb", "A"
        )
        self.assertEqual(info["native_seq"][:5], "NLYIQ")
        self.assertEqual(len(info["resseqs"]), 20)


if __name__ == "__main__":
    unittest.main()
