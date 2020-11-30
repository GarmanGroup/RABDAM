
# python -m unittest tests/test_atom_filtering.py

import unittest
from rabdam.Subroutines.parsePDB import atom, b_damage_atom_list
from rabdam.Subroutines.PDBCUR import clean_atom_rec

class TestClass(unittest.TestCase):

    def test_clean_atom_records(self):
        """
        Checks that ATOM/HETATM records are filtered appropriately to remove
        hydrogen atoms, 0 occupancy atoms, retain only a single conformer
        per-residue, and check that disulfide bonds have not been refined to
        sub-1 occupancy
        """

    def test_create_list_of_atoms_for_bdamage_analysis(self):
        """
        Checks that input list of atom objects is correctly filtered
        according to user-specified program inputs
        """
