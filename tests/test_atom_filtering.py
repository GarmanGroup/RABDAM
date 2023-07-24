
# RABDAM
# Copyright (C) 2023 Garman Group, University of Oxford

# This file is part of RABDAM.

# RABDAM is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation, either version 3 of
# the License, or (at your option) any later version.

# RABDAM is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General
# Public License along with this program.  If not, see
# <http://www.gnu.org/licenses/>.


# python -m unittest tests/test_atom_filtering.py

import unittest
from rabdam.Subroutines.parsePDB import atom, b_damage_atom_list
from rabdam.Subroutines.PDBCUR import clean_atom_rec

class TestClass(unittest.TestCase):

    def make_atoms_list(self):
        """
        """

        atoms_list = []
        # Taken from 2BN1 and 2QTZ
        new_atom_1 = atom()
        new_atom_1.lineID = 'ATOM'
        new_atom_1.atomNum = 1
        new_atom_1.atomType = 'CA'
        new_atom_1.conformer = ''
        new_atom_1.resiType = 'GLY'
        new_atom_1.chainID = 'A'
        new_atom_1.resiNum = 1
        new_atom_1.insCode = ''
        new_atom_1.xyzCoords = [[13.602], [45.768], [30.728]]
        new_atom_1.occupancy = 1.0
        new_atom_1.bFactor = 19.08
        new_atom_1.element = 'C'
        new_atom_1.charge = ''
        new_atom_1.origResiNum = ''
        new_atom_1.origResiType = 'GLY'
        new_atom_1.origChainID = 'A'
        new_atom_1.origAtomType = 'CA'
        new_atom_1.pd = 19
        new_atom_1.avrg_bf = 26.2
        new_atom_1.bd = 1.3
        new_atom_1.protein = True
        new_atom_1.na = False
        new_atom_1.chain_len = 100
        atoms_list.append(new_atom_1)

        new_atom_2 = atom()
        new_atom_2.lineID = 'ATOM'
        new_atom_2.atomNum = 2
        new_atom_2.atomType = 'CA'
        new_atom_2.conformer = ''
        new_atom_2.resiType = 'GLY'
        new_atom_2.chainID = 'A'
        new_atom_2.resiNum = 2
        new_atom_2.insCode = ''
        new_atom_2.xyzCoords = [[13.128], [46.298], [34.342]]
        new_atom_2.occupancy = 1.0
        new_atom_2.bFactor = 33.1
        new_atom_2.element = 'C'
        new_atom_2.charge = ''
        new_atom_2.origResiNum = ''
        new_atom_2.origResiType = 'GLY'
        new_atom_2.origChainID = 'A'
        new_atom_2.origAtomType = 'CA'
        new_atom_2.pd = 30
        new_atom_2.avrg_bf = 25.6
        new_atom_2.bd = 0.97
        new_atom_2.protein = True
        new_atom_2.na = False
        new_atom_2.chain_len = 100
        atoms_list.append(new_atom_2)

        new_atom_3 = atom()
        new_atom_3.lineID = 'HETATM'
        new_atom_3.atomNum = 3
        new_atom_3.atomType = 'PA'
        new_atom_3.conformer = 'X'
        new_atom_3.resiType = 'FAD'
        new_atom_3.chainID = 'B'
        new_atom_3.resiNum = 2
        new_atom_3.insCode = ''
        new_atom_3.xyzCoords = [[29.643], [29.7], [51.402]]
        new_atom_3.occupancy = 1.0
        new_atom_3.bFactor = 20.52
        new_atom_3.element = 'P'
        new_atom_3.charge = ''
        new_atom_3.origResiNum = ''
        new_atom_3.origResiType = 'FAD'
        new_atom_3.origChainID = 'B'
        new_atom_3.origAtomType = 'PA'
        new_atom_3.pd = 26
        new_atom_3.avrg_bf = 20.5
        new_atom_3.bd = 2.0
        new_atom_3.protein = False
        new_atom_3.na = False
        new_atom_3.chain_len = 100
        atoms_list.append(new_atom_3)

        return atoms_list
    
    def test_clean_atom_records(self):
        """
        Checks that ATOM/HETATM records are filtered appropriately to remove
        hydrogen atoms, 0 occupancy atoms, retain only a single conformer
        per-residue, and check that disulfide bonds have not been refined to
        sub-1 occupancy
        """

        import copy
        import os
        import shutil

        if not os.path.isdir('tests/temp_files/'):
            os.mkdir('tests/temp_files/')

        atoms_list_1 = self.make_atoms_list()

        # Should successfully parse all atoms
        exp_atoms_list_1 = copy.deepcopy(atoms_list_1)
        exp_file_content_1 = [
            '',
            'ATOM      1  CA  GLY A   1      13.602  45.768  30.728  1.00 19.08           C  ',
            'ATOM      2  CA  GLY A   2      13.128  46.298  34.342  1.00 33.10           C  ',
            'TER                                                                             ',
            'HETATM    3  PA XFAD B   2      29.643  29.700  51.402  1.00 20.52           P  ',
            'END                                                                             ',
            ''
        ]
        exit_1, pause_1, act_atoms_list_1, file_1 = clean_atom_rec(
            atoms_list_1, '', 'tests/temp_files/test'
        )
        self.assertFalse(exit_1)
        self.assertFalse(pause_1)
        self.assertEqual(exp_atoms_list_1, act_atoms_list_1)
        with open('tests/temp_files/test_asymmetric_unit.pdb', 'r') as f:
            act_file_content_1 = f.read().split('\n')
            print(act_file_content_1)
            self.assertEqual(exp_file_content_1, act_file_content_1)

        # Hydrogen atoms - should remove atom 1 when set to hydrogen
        atoms_list_2 = copy.deepcopy(atoms_list_1)
        atoms_list_2[0].element = 'H'
        exp_atoms_list_2 = copy.deepcopy(atoms_list_2)[1:]
        exp_file_content_2 = [
            '',
            'ATOM      2  CA  GLY A   2      13.128  46.298  34.342  1.00 33.10           C  ',
            'TER                                                                             ',
            'HETATM    3  PA XFAD B   2      29.643  29.700  51.402  1.00 20.52           P  ',
            'END                                                                             ',
            ''
        ]
        exit_2, pause_2, act_atoms_list_2, file_2 = clean_atom_rec(
            atoms_list_2, '', 'tests/temp_files/test'
        )
        self.assertFalse(exit_2)
        self.assertFalse(pause_2)
        self.assertEqual(exp_atoms_list_2, act_atoms_list_2)
        with open('tests/temp_files/test_asymmetric_unit.pdb', 'r') as f:
            act_file_content_2 = f.read().split('\n')
            self.assertEqual(exp_file_content_2, act_file_content_2)

        # 0 occupancy - should remove atom 2 when occupancy set to 0
        atoms_list_3 = copy.deepcopy(atoms_list_1)
        atoms_list_3[1].occupancy = 0
        exp_atoms_list_3 = [copy.deepcopy(atoms_list_3)[0], copy.deepcopy(atoms_list_3)[2]]
        exp_file_content_3 = [
            '',
            'ATOM      1  CA  GLY A   1      13.602  45.768  30.728  1.00 19.08           C  ',
            'TER                                                                             ',
            'HETATM    3  PA XFAD B   2      29.643  29.700  51.402  1.00 20.52           P  ',
            'END                                                                             ',
            ''
        ]
        exit_3, pause_3, act_atoms_list_3, file_3 = clean_atom_rec(
            atoms_list_3, '', 'tests/temp_files/test'
        )
        self.assertFalse(exit_3)
        self.assertFalse(pause_3)
        self.assertEqual(exp_atoms_list_3, act_atoms_list_3)
        with open('tests/temp_files/test_asymmetric_unit.pdb', 'r') as f:
            act_file_content_3 = f.read().split('\n')
            self.assertEqual(exp_file_content_3, act_file_content_3)

        # Negative Bfactor - should remove atom 3 when its Bfactor is set to a
        # value < 0
        atoms_list_4 = copy.deepcopy(atoms_list_1)
        atoms_list_4[2].bFactor = -0.1
        exp_atoms_list_4 = copy.deepcopy(atoms_list_4)[:2]
        exp_file_content_4 = [
            '',
            'ATOM      1  CA  GLY A   1      13.602  45.768  30.728  1.00 19.08           C  ',
            'ATOM      2  CA  GLY A   2      13.128  46.298  34.342  1.00 33.10           C  ',
            'TER                                                                             ',
            'END                                                                             ',
            ''
        ]
        exit_4, pause_4, act_atoms_list_4, file_4 = clean_atom_rec(
            atoms_list_4, '', 'tests/temp_files/test'
        )
        self.assertFalse(exit_4)
        self.assertFalse(pause_4)
        self.assertEqual(exp_atoms_list_4, act_atoms_list_4)
        with open('tests/temp_files/test_asymmetric_unit.pdb', 'r') as f:
            act_file_content_4 = f.read().split('\n')
            self.assertEqual(exp_file_content_4, act_file_content_4)

        # Single sub-1 occupancy - should pause the program (whilst retaining
        # all atoms)
        atoms_list_5 = copy.deepcopy(atoms_list_1)
        atoms_list_5[2].occupancy = 0.7
        exp_atoms_list_5 = copy.deepcopy(atoms_list_5)
        exp_file_content_5 = [
            '',
            'ATOM      1  CA  GLY A   1      13.602  45.768  30.728  1.00 19.08           C  ',
            'ATOM      2  CA  GLY A   2      13.128  46.298  34.342  1.00 33.10           C  ',
            'TER                                                                             ',
            'HETATM    3  PA XFAD B   2      29.643  29.700  51.402  0.70 20.52           P  ',
            'END                                                                             ',
            ''
        ]
        exit_5, pause_5, act_atoms_list_5, file_5 = clean_atom_rec(
            atoms_list_5, '', 'tests/temp_files/test'
        )
        self.assertFalse(exit_5)
        self.assertTrue(pause_5)
        self.assertEqual(exp_atoms_list_5, act_atoms_list_5)
        with open('tests/temp_files/test_asymmetric_unit.pdb', 'r') as f:
            act_file_content_5 = f.read().split('\n')
            self.assertEqual(exp_file_content_5, act_file_content_5)

        # Alternate conformers that sum to 1 - should retain the highest
        #  occupancy conformer
        atoms_list_6 = copy.deepcopy(atoms_list_1)
        atoms_list_6[0].conformer = 'A'
        atoms_list_6[0].occupancy = 0.25
        atoms_list_6[1].resiNum = 1
        atoms_list_6[1].conformer = 'B'
        atoms_list_6[1].occupancy = 0.75
        exp_atoms_list_6 = copy.deepcopy(atoms_list_6)[1:]
        exp_file_content_6 = [
            '',
            'ATOM      2  CA BGLY A   1      13.128  46.298  34.342  0.75 33.10           C  ',
            'TER                                                                             ',
            'HETATM    3  PA XFAD B   2      29.643  29.700  51.402  1.00 20.52           P  ',
            'END                                                                             ',
            ''
        ]
        exit_6, pause_6, act_atoms_list_6, file_6 = clean_atom_rec(
            atoms_list_6, '', 'tests/temp_files/test'
        )
        self.assertFalse(exit_6)
        self.assertFalse(pause_6)
        self.assertEqual(exp_atoms_list_6, act_atoms_list_6)
        with open('tests/temp_files/test_asymmetric_unit.pdb', 'r') as f:
            act_file_content_6 = f.read().split('\n')
            self.assertEqual(exp_file_content_6, act_file_content_6)

        # Alternate conformers that don't sum to 1 - should pause the program
        # (whilst retaining only the highest occupancy conformer)
        atoms_list_7 = copy.deepcopy(atoms_list_1)
        atoms_list_7[0].conformer = 'A'
        atoms_list_7[0].occupancy = 0.25
        atoms_list_7[1].resiNum = 1
        atoms_list_7[1].conformer = 'B'
        atoms_list_7[1].occupancy = 0.74
        exp_atoms_list_7 = copy.deepcopy(atoms_list_7)[1:]
        exp_file_content_7 = [
            '',
            'ATOM      2  CA BGLY A   1      13.128  46.298  34.342  0.74 33.10           C  ',
            'TER                                                                             ',
            'HETATM    3  PA XFAD B   2      29.643  29.700  51.402  1.00 20.52           P  ',
            'END                                                                             ',
            ''
        ]
        exit_7, pause_7, act_atoms_list_7, file_7 = clean_atom_rec(
            atoms_list_7, '', 'tests/temp_files/test'
        )
        self.assertFalse(exit_7)
        self.assertTrue(pause_7)
        self.assertEqual(exp_atoms_list_7, act_atoms_list_7)
        with open('tests/temp_files/test_asymmetric_unit.pdb', 'r') as f:
            act_file_content_7 = f.read().split('\n')
            self.assertEqual(exp_file_content_7, act_file_content_7)

        # No atoms left after filtering - should set exit to True
        atoms_list_8 = copy.deepcopy(atoms_list_1)
        atoms_list_8[0].element = 'H'
        atoms_list_8[1].occupancy = 1.4
        atoms_list_8[2].bFactor = 0
        exp_atoms_list_8 = []
        exp_file_content_8 = [
            '',
            'END                                                                             ',
            ''
        ]
        exit_8, pause_8, act_atoms_list_8, file_8 = clean_atom_rec(
            atoms_list_8, '', 'tests/temp_files/test'
        )
        self.assertTrue(exit_8)
        self.assertFalse(pause_8)
        self.assertEqual(exp_atoms_list_8, act_atoms_list_8)
        with open('tests/temp_files/test_asymmetric_unit.pdb', 'r') as f:
            act_file_content_8 = f.read().split('\n')
            self.assertEqual(exp_file_content_8, act_file_content_8)

        shutil.rmtree('tests/temp_files/')

    def test_create_list_of_atoms_for_bdamage_analysis(self):
        """
        Checks that input list of atom objects is correctly filtered
        according to user-specified program inputs
        """

        import copy

        atoms_list_1 = self.make_atoms_list()

        # Remove/keep HETATM
        # Set to remove HETATM - should remove atom 3
        act_atoms_list_1 = b_damage_atom_list(
            atoms_list_1, False, 'protein', [], []
        )
        exp_atoms_list_1 = copy.deepcopy(atoms_list_1)[:2]
        self.assertEqual(exp_atoms_list_1, act_atoms_list_1)

        # Set to remove HETATM - should keep all three atoms since atom is now
        # defined as "protein"
        atoms_list_2 = copy.deepcopy(atoms_list_1)
        atoms_list_2[2].protein = True
        act_atoms_list_2 = b_damage_atom_list(
            atoms_list_2, False, 'protein', [], []
        )
        exp_atoms_list_2 = copy.deepcopy(atoms_list_2)
        self.assertEqual(exp_atoms_list_2, act_atoms_list_2)

        # Set to remove HETATM - should remove atom 3 since atom is now
        # defined as "na", and program is directed to remove non-protein atoms
        atoms_list_3 = copy.deepcopy(atoms_list_1)
        atoms_list_3[2].na = True
        act_atoms_list_3 = b_damage_atom_list(
            atoms_list_3, False, 'protein', [], []
        )
        exp_atoms_list_3 = copy.deepcopy(atoms_list_3)[:2]
        self.assertEqual(exp_atoms_list_3, act_atoms_list_3)

        # Set to remove HETATM - should keep atom 3 since atom 3 is now
        # defined as "na", and program is directed to keep nucleic acid atoms
        atoms_list_4 = copy.deepcopy(atoms_list_1)
        atoms_list_4[2].na = True
        act_atoms_list_4 = b_damage_atom_list(
            atoms_list_4, False, 'na', [], []
        )
        exp_atoms_list_4 = copy.deepcopy(atoms_list_4)[2:]
        self.assertEqual(exp_atoms_list_4, act_atoms_list_4)

        # Set to remove HETATM - should remove all atoms since are all
        # defined as "protein", and program is directed to remove non-NA atoms
        atoms_list_5 = copy.deepcopy(atoms_list_1)
        atoms_list_5[2].protein = True
        act_atoms_list_5 = b_damage_atom_list(
            atoms_list_5, False, 'na', [], []
        )
        exp_atoms_list_5 = []
        self.assertEqual(exp_atoms_list_5, act_atoms_list_5)

        # Set to keep HETATM - should retain all atoms
        atoms_list_6 = copy.deepcopy(atoms_list_1)
        act_atoms_list_6 = b_damage_atom_list(
            atoms_list_6, True, 'protein', [], []
        )
        exp_atoms_list_6 = copy.deepcopy(atoms_list_6)
        self.assertEqual(exp_atoms_list_6, act_atoms_list_6)

        # Set to remove HETATM, and keep all protein and NA atoms.
        atoms_list_7 = copy.deepcopy(atoms_list_1)
        atoms_list_7[1].na = True
        atoms_list_7[1].protein = False
        act_atoms_list_7 = b_damage_atom_list(
            atoms_list_7, False, 'proteinna', [], []
        )
        exp_atoms_list_7 = copy.deepcopy(atoms_list_7)[:2]
        self.assertEqual(exp_atoms_list_7, act_atoms_list_7)

        # Remove atoms -should remove atoms 1 and 3
        atoms_list_8 = copy.deepcopy(atoms_list_1)
        act_atoms_list_8 = b_damage_atom_list(
            atoms_list_8, True, 'proteinna', [], ['1', '3']
        )
        exp_atoms_list_8 = copy.deepcopy(atoms_list_8)[1:2]
        self.assertEqual(exp_atoms_list_8, act_atoms_list_8)

        # Should remove all atoms in GLY residues
        atoms_list_9 = copy.deepcopy(atoms_list_1)
        act_atoms_list_9 = b_damage_atom_list(
            atoms_list_9, True, 'proteinna', [], ['GLY']
        )
        exp_atoms_list_9 = copy.deepcopy(atoms_list_9)[2:]
        self.assertEqual(exp_atoms_list_9, act_atoms_list_9)

        # Should remove atoms 1 and 3, then add back in any atoms in GLY residues
        atoms_list_10 = copy.deepcopy(atoms_list_1)
        act_atoms_list_10 = b_damage_atom_list(
            atoms_list_10, True, 'proteinna', ['GLY'], ['1', '3']
        )
        exp_atoms_list_10 = copy.deepcopy(atoms_list_10)[:2]
        self.assertEqual(exp_atoms_list_10, act_atoms_list_10)
