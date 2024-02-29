
# RABDAM
# Copyright (C) 2024 Garman Group, University of Oxford

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
from tests.make_atom_rec import gen_atom_objs_list

class TestClass(unittest.TestCase):
    
    def test_clean_atom_rec(self):
        """
        Checks that ATOM/HETATM records are filtered appropriately to remove
        hydrogen atoms, 0 occupancy atoms, retain only a single conformer
        per-residue, and check that Glu and Asp residues have not been refined
        to sub-1 occupancy (across all conformers)
        """

        import copy
        import os
        import shutil

        if not os.path.isdir('tests/temp_files/'):
            os.mkdir('tests/temp_files/')

        atoms_list_1 = gen_atom_objs_list()
        cif_header = [
            'data_tests/temp_files/test_asymmetric_unit',
            '#',
            'loop_',
            '_atom_site.group_PDB',
            '_atom_site.id',
            '_atom_site.type_symbol',
            '_atom_site.label_atom_id',
            '_atom_site.label_alt_id',
            '_atom_site.label_comp_id',
            '_atom_site.label_asym_id',
            '_atom_site.label_entity_id',
            '_atom_site.label_seq_id',
            '_atom_site.pdbx_PDB_ins_code',
            '_atom_site.Cartn_x',
            '_atom_site.Cartn_y',
            '_atom_site.Cartn_z',
            '_atom_site.occupancy',
            '_atom_site.B_iso_or_equiv',
            '_atom_site.pdbx_formal_charge',
            '_atom_site.auth_seq_id',
            '_atom_site.auth_comp_id',
            '_atom_site.auth_asym_id',
            '_atom_site.auth_atom_id',
            '_atom_site.pdbx_PDB_model_num'
        ]

        # Should successfully parse all atoms
        exp_atoms_list_1 = copy.deepcopy(atoms_list_1)
        exp_cif_lines_1 = cif_header + [
            'ATOM   2    C CA . GLY A 1 1 ? 13.602 45.768 30.728 1.0 19.08 ? 1    GLY A CA 1.0',
            'HETATM 438  O O  . HOH C 2 . ? 13.128 46.298 34.342 1.0 33.1  ? 2001 HOH B O  1.1',
            'HETATM 3546 P PA X FAD B 3 . ? 29.643 29.7   51.402 1.0 20.52 ? 700  FAD B PA 12 ',
            '#',
            ''
        ]
        (
            exit_1, pause_1, act_atoms_list_1, clean_au_file_1,
            sub_1_asp_glu_occ_1
        ) = clean_atom_rec(atoms_list_1, 'tests/temp_files/test', False)
        self.assertFalse(exit_1)
        self.assertFalse(pause_1)
        self.assertEqual(exp_atoms_list_1, act_atoms_list_1)
        with open('{}.cif'.format(clean_au_file_1), 'r') as f:
            act_cif_lines_1 = f.read().split('\n')
            self.assertEqual(exp_cif_lines_1, act_cif_lines_1)
        self.assertFalse(sub_1_asp_glu_occ_1)

        # Check removes atom 1 when it is set to a hydrogen atom
        atoms_list_2 = copy.deepcopy(atoms_list_1)
        atoms_list_2[0].element = 'H'
        exp_atoms_list_2 = copy.deepcopy(atoms_list_2)[1:]
        exp_cif_lines_2 = cif_header + [
            'HETATM 438  O O  . HOH C 2 ? ? 13.128 46.298 34.342 1.0 33.1  ? 2001 HOH B O  1.1',
            'HETATM 3546 P PA X FAD B 3 ? ? 29.643 29.7   51.402 1.0 20.52 ? 700  FAD B PA 12 ',
            '#',
            ''
        ]
        (
            exit_2, pause_2, act_atoms_list_2, clean_au_file_2,
            sub_1_asp_glu_occ_2
        ) = clean_atom_rec(atoms_list_2, 'tests/temp_files/test', False)
        self.assertFalse(exit_2)
        self.assertFalse(pause_2)
        self.assertEqual(exp_atoms_list_2, act_atoms_list_2)
        with open('{}.cif'.format(clean_au_file_2), 'r') as f:
            act_cif_lines_2 = f.read().split('\n')
            self.assertEqual(exp_cif_lines_2, act_cif_lines_2)
        self.assertFalse(sub_1_asp_glu_occ_2)

        # Check removes atom 2 when occupancy is set to 0
        atoms_list_3 = copy.deepcopy(atoms_list_1)
        atoms_list_3[1].occupancy = 0.0
        exp_atoms_list_3 = [copy.deepcopy(atoms_list_3)[0], copy.deepcopy(atoms_list_3)[2]]
        exp_cif_lines_3 = cif_header + [
            'ATOM   2    C CA . GLY A 1 1 ? 13.602 45.768 30.728 1.0 19.08 ? 1   GLY A CA 1.0',
            'HETATM 3546 P PA X FAD B 3 . ? 29.643 29.7   51.402 1.0 20.52 ? 700 FAD B PA 12 ',
            '#',
            ''
        ]
        (
            exit_3, pause_3, act_atoms_list_3, clean_au_file_3,
            sub_1_asp_glu_occ_3
        ) = clean_atom_rec(atoms_list_3, 'tests/temp_files/test', False)
        self.assertFalse(exit_3)
        self.assertFalse(pause_3)
        self.assertEqual(exp_atoms_list_3, act_atoms_list_3)
        with open('{}.cif'.format(clean_au_file_3), 'r') as f:
            act_cif_lines_3 = f.read().split('\n')
            self.assertEqual(exp_cif_lines_3, act_cif_lines_3)
        self.assertFalse(sub_1_asp_glu_occ_3)

        # Check removes atom 3 when Bfactor <= 0
        atoms_list_4 = copy.deepcopy(atoms_list_1)
        atoms_list_4[2].bFactor = 0
        exp_atoms_list_4 = copy.deepcopy(atoms_list_4)[:2]
        exp_cif_lines_4 = cif_header + [
            'ATOM   2   C CA ? GLY A 1 1 ? 13.602 45.768 30.728 1.0 19.08 ? 1    GLY A CA 1.0',
            'HETATM 438 O O  ? HOH C 2 . ? 13.128 46.298 34.342 1.0 33.1  ? 2001 HOH B O  1.1',
            '#',
            ''
        ]
        (
            exit_4, pause_4, act_atoms_list_4, clean_au_file_4,
            sub_1_asp_glu_occ_4
        ) = clean_atom_rec(atoms_list_4, 'tests/temp_files/test', False)
        self.assertFalse(exit_4)
        self.assertFalse(pause_4)
        self.assertEqual(exp_atoms_list_4, act_atoms_list_4)
        with open('{}.cif'.format(clean_au_file_4), 'r') as f:
            act_cif_lines_4 = f.read().split('\n')
            self.assertEqual(exp_cif_lines_4, act_cif_lines_4)
        self.assertFalse(sub_1_asp_glu_occ_4)

        # Check pauses program if non-ASP/GLU atom has sub-1 occupancy
        atoms_list_5 = copy.deepcopy(atoms_list_1)
        atoms_list_5[2].occupancy = 0.7
        exp_atoms_list_5 = copy.deepcopy(atoms_list_5)
        exp_cif_lines_5 = cif_header + [
            'ATOM   2    C CA . GLY A 1 1 ? 13.602 45.768 30.728 1.0 19.08 ? 1    GLY A CA 1.0',
            'HETATM 438  O O  . HOH C 2 . ? 13.128 46.298 34.342 1.0 33.1  ? 2001 HOH B O  1.1',
            'HETATM 3546 P PA X FAD B 3 . ? 29.643 29.7   51.402 0.7 20.52 ? 700  FAD B PA 12 ',
            '#',
            ''
        ]
        (
            exit_5, pause_5, act_atoms_list_5, clean_au_file_5,
            sub_1_asp_glu_occ_5
        ) = clean_atom_rec(atoms_list_5, 'tests/temp_files/test', False)
        self.assertFalse(exit_5)
        self.assertTrue(pause_5)
        self.assertEqual(exp_atoms_list_5, act_atoms_list_5)
        with open('{}.cif'.format(clean_au_file_5), 'r') as f:
            act_cif_lines_5 = f.read().split('\n')
            self.assertEqual(exp_cif_lines_5, act_cif_lines_5)
        self.assertFalse(sub_1_asp_glu_occ_5)

        # Check pauses program and sets sub_1_asp_glu_occ to True if
        # non-ASP/GLU atom has sub-1 occupancy
        atoms_list_6 = copy.deepcopy(atoms_list_1)
        atoms_list_6[2].resiType = 'GLU'
        atoms_list_6[2].origResiType = 'GLU'
        atoms_list_6[2].occupancy = 0.7
        exp_atoms_list_6 = copy.deepcopy(atoms_list_6)
        exp_cif_lines_6 = cif_header + [
            'ATOM   2    C CA . GLY A 1 1 ? 13.602 45.768 30.728 1.0 19.08 ? 1    GLY A CA 1.0',
            'HETATM 438  O O  . HOH C 2 . ? 13.128 46.298 34.342 1.0 33.1  ? 2001 HOH B O  1.1',
            'HETATM 3546 P PA X GLU B 3 . ? 29.643 29.7   51.402 0.7 20.52 ? 700  GLU B PA 12 ',
            '#',
            ''
        ]
        (
            exit_6, pause_6, act_atoms_list_6, clean_au_file_6,
            sub_1_asp_glu_occ_6
        ) = clean_atom_rec(atoms_list_6, 'tests/temp_files/test', False)
        self.assertFalse(exit_6)
        self.assertTrue(pause_6)
        self.assertEqual(exp_atoms_list_6, act_atoms_list_6)
        with open('{}.cif'.format(clean_au_file_6), 'r') as f:
            act_cif_lines_6 = f.read().split('\n')
            self.assertEqual(exp_cif_lines_6, act_cif_lines_6)
        self.assertTrue(sub_1_asp_glu_occ_6)

        # Check retains highest occupancy conformer when occupancies sum to 1
        atoms_list_7 = copy.deepcopy(atoms_list_1)
        atoms_list_7[0].conformer = 'A'
        atoms_list_7[0].occupancy = 0.25

        atoms_list_7[1].chainID = 'A'
        atoms_list_7[1].resiNum = 1
        atoms_list_7[1].resiType = 'GLY'
        atoms_list_7[1].atomType = 'CA'
        atoms_list_7[1].conformer = 'B'
        atoms_list_7[1].occupancy = 0.75

        exp_atoms_list_7 = copy.deepcopy(atoms_list_7[1:])
        exp_cif_lines_7 = cif_header + [
            'HETATM 438  O O  B HOH C 2 ? ? 13.128 46.298 34.342 0.75 33.1  ? 1   GLY A CA 1.1',
            'HETATM 3546 P PA X FAD B 3 ? ? 29.643 29.7   51.402 1.0  20.52 ? 700 FAD B PA 12 ',
            '#',
            ''
        ]
        (
            exit_7, pause_7, act_atoms_list_7, clean_au_file_7,
            sub_1_asp_glu_occ_7
        ) = clean_atom_rec(atoms_list_7, 'tests/temp_files/test', False)
        self.assertFalse(exit_7)
        self.assertFalse(pause_7)
        self.assertEqual(exp_atoms_list_7, act_atoms_list_7)
        with open('{}.cif'.format(clean_au_file_7), 'r') as f:
            act_cif_lines_7 = f.read().split('\n')
            self.assertEqual(exp_cif_lines_7, act_cif_lines_7)
        self.assertFalse(sub_1_asp_glu_occ_7)

        # Check pauses program when occupancies of alternate conformers do not
        # sum to 1
        atoms_list_8 = copy.deepcopy(atoms_list_1)
        atoms_list_8[0].conformer = 'A'
        atoms_list_8[0].occupancy = 0.25
        atoms_list_8[0].resiType = 'DAS'

        atoms_list_8[1].chainID = 'A'
        atoms_list_8[1].resiNum = 1
        atoms_list_8[1].resiType = 'DAS'
        atoms_list_8[1].atomType = 'CA'
        atoms_list_8[1].conformer = 'B'
        atoms_list_8[1].occupancy = 0.74

        exp_atoms_list_8 = copy.deepcopy(atoms_list_8[1:])
        exp_cif_lines_8 = cif_header + [
            'HETATM 438  O O  B HOH C 2 ? ? 13.128 46.298 34.342 0.74 33.1  ? 1   DAS A CA 1.1',
            'HETATM 3546 P PA X FAD B 3 ? ? 29.643 29.7   51.402 1.0  20.52 ? 700 FAD B PA 12 ',
            '#',
            ''
        ]
        (
            exit_8, pause_8, act_atoms_list_8, clean_au_file_8,
            sub_1_asp_glu_occ_8
        ) = clean_atom_rec(atoms_list_8, 'tests/temp_files/test', False)
        self.assertFalse(exit_8)
        self.assertTrue(pause_8)
        self.assertEqual(exp_atoms_list_8, act_atoms_list_8)
        with open('{}.cif'.format(clean_au_file_8), 'r') as f:
            act_cif_lines_8 = f.read().split('\n')
            self.assertEqual(exp_cif_lines_8, act_cif_lines_8)
        self.assertTrue(sub_1_asp_glu_occ_8)

        # Check sets exit to True when no atoms remain after filtering
        atoms_list_9 = copy.deepcopy(atoms_list_1)
        atoms_list_9[0].occupancy = 1.1
        atoms_list_9[1].bFactor = -20
        atoms_list_9[2].element = 'H'

        exp_atoms_list_9 = []
        exp_cif_lines_9 = cif_header + [
            '#',
            ''
        ]
        (
            exit_9, pause_9, act_atoms_list_9, clean_au_file_9,
            sub_1_asp_glu_occ_9
        ) = clean_atom_rec(atoms_list_9, 'tests/temp_files/test', False)
        self.assertTrue(exit_9)
        self.assertFalse(pause_9)
        self.assertEqual(exp_atoms_list_9, act_atoms_list_9)
        with open('{}.cif'.format(clean_au_file_9), 'r') as f:
            act_cif_lines_9 = f.read().split('\n')
            self.assertEqual(exp_cif_lines_9, act_cif_lines_9)
        self.assertFalse(sub_1_asp_glu_occ_9)

        shutil.rmtree('tests/temp_files/')

    def test_b_damage_atom_list(self):
        """
        Checks that input list of atom objects is correctly filtered according
        to user-specified program inputs
        """

        import copy

        atoms_list = gen_atom_objs_list()

        """
        'ATOM   2    C CA . GLY A 1 1 ? 13.602 45.768 30.728 1.0 19.08 ? 1    GLY A CA 1.0',
        'HETATM 438  O O  . HOH C 2 . ? 13.128 46.298 34.342 1.0 33.1  ? 2001 HOH B O  1.1',
        'HETATM 3546 P PA X FAD B 3 . ? 29.643 29.7   51.402 1.0 20.52 ? 700  FAD B PA 12 ',
        """

        # Check discards atoms where protein=False and na=False when HETATM set
        # to False
        atoms_list_1 = copy.deepcopy(atoms_list)
        act_atoms_list_1 = b_damage_atom_list(
            clean_au_list=atoms_list_1, HETATM=False, protOrNA='proteinna',
            addAtoms=[], removeAtoms=[]
        )
        exp_atoms_list_1 = copy.deepcopy(atoms_list_1)[1:]
        self.assertEqual(exp_atoms_list_1, act_atoms_list_1)

        # Check keeps all atoms where protein=False and na=False when HETATM set
        # to True
        atoms_list_2 = copy.deepcopy(atoms_list)
        act_atoms_list_2 = b_damage_atom_list(
            clean_au_list=atoms_list_2, HETATM=True, protOrNA='proteinna',
            addAtoms=[], removeAtoms=[]
        )
        exp_atoms_list_2 = copy.deepcopy(atoms_list_2)
        self.assertEqual(exp_atoms_list_2, act_atoms_list_2)

        # Check removes nucleic acids if protOrNA is set to 'protein'
        atoms_list_3 = copy.deepcopy(atoms_list)
        act_atoms_list_3 = b_damage_atom_list(
            clean_au_list=atoms_list_3, HETATM=True, protOrNA='protein',
            addAtoms=[], removeAtoms=[]
        )
        exp_atoms_list_3 = copy.deepcopy(atoms_list_3[:2])
        self.assertEqual(exp_atoms_list_3, act_atoms_list_3)

        # Check removes protein if protOrNA is set to 'na'
        atoms_list_4 = copy.deepcopy(atoms_list)
        atoms_list_4[0].protein = True
        act_atoms_list_4 = b_damage_atom_list(
            clean_au_list=atoms_list_4, HETATM=True, protOrNA='protein',
            addAtoms=[], removeAtoms=[]
        )
        exp_atoms_list_4 = copy.deepcopy(atoms_list_4[:2])
        self.assertEqual(exp_atoms_list_4, act_atoms_list_4)

        # Check removes atoms 1 and 3 if remove_atoms = ['2', '3546']
        atoms_list_5 = copy.deepcopy(atoms_list)
        act_atoms_list_5 = b_damage_atom_list(
            clean_au_list=atoms_list_5, HETATM=True, protOrNA='proteinna',
            addAtoms=[], removeAtoms=['2', '3546']
        )
        exp_atoms_list_5 = copy.deepcopy(atoms_list_5[1:2])
        self.assertEqual(exp_atoms_list_5, act_atoms_list_5)

        # Check removes atom 1 if remove_atoms = ['GLY']
        atoms_list_6 = copy.deepcopy(atoms_list)
        act_atoms_list_6 = b_damage_atom_list(
            clean_au_list=atoms_list_6, HETATM=True, protOrNA='proteinna',
            addAtoms=[], removeAtoms=['GLY']
        )
        exp_atoms_list_6 = copy.deepcopy(atoms_list_6[1:])
        self.assertEqual(exp_atoms_list_6, act_atoms_list_6)

        # Check adds atom 1 back in but doesn't duplicate atoms 1 or 2 if
        # add_atoms = ['2', '438', 'GLY']
        atoms_list_7 = copy.deepcopy(atoms_list)
        atoms_list_7[1].protein = False
        act_atoms_list_7 = b_damage_atom_list(
            clean_au_list=atoms_list_7, HETATM=False, protOrNA='na',
            addAtoms=[], removeAtoms=[]
        )
        exp_atoms_list_7 = copy.deepcopy(atoms_list_7[1:])
        self.assertEqual(exp_atoms_list_7, act_atoms_list_7)

        atoms_list_8 = copy.deepcopy(atoms_list)
        atoms_list_8[1].protein = False
        act_atoms_list_8 = b_damage_atom_list(
            clean_au_list=atoms_list_8, HETATM=False, protOrNA='na',
            addAtoms=['2', '438', 'GLY'], removeAtoms=[]
        )
        exp_atoms_list_8 = copy.deepcopy(atoms_list_8)
        self.assertEqual(exp_atoms_list_8, act_atoms_list_8)

