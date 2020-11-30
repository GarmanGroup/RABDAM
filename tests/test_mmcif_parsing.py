
# python -m unittest tests/test_mmcif_parsing.py

import unittest
from rabdam.Subroutines.parsePDB import atom, download_mmcif
from rabdam.Subroutines.PDBCUR import (
    make_cryst1_line_from_mmcif, find_disulfides_from_mmcif,
    parse_seqres_from_mmcif, parse_atom_rec_from_mmcif
)
from rabdam.Subroutines.output import generate_output_files

class TestClass(unittest.TestCase):

    def test_download_mmcif(self):
        """
        Checks that downloads mmCIF files correctly, and raises the appropriate
        errors e.g. if the file does not exist
        """

        import os
        import shutil

        exp_results = {'2BN1': False,
                       '1ABC': True}
        for pdb, exp_exit in exp_results.items():
            if os.path.isdir('tests/temp_files/'):
                shutil.rmtree('tests/temp_files/')
            act_exit = download_mmcif(
                pdb, 'tests/temp_files/', 'tests/temp_files/test.cif'
            )
            self.assertEqual(exp_exit, act_exit)

    def test_make_cryst1_line_from_mmcif(self):
        """
        Checks that CRYST1 line (found in PDB format files, and required for
        running PDBCUR to generate the unit cell) is generated correctly from an
        mmCIF format input file
        """

        import copy

        # Taken from 2BN1
        space_group_info_1 = [
            '_cell.entry_id           2BN1',
            '_cell.length_a           77.900',
            '_cell.length_b           77.900',
            '_cell.length_c           77.900',
            '_cell.angle_alpha        90.00',
            '_cell.angle_beta         90.00',
            '_cell.angle_gamma        90.00',
            '_cell.Z_PDB              24',
            '_cell.pdbx_unique_axis   ?',
            '_symmetry.entry_id                         2BN1',
            '_symmetry.space_group_name_H-M             \'I 21 3\'',
            '_symmetry.pdbx_full_space_group_name_H-M   ?',
            '_symmetry.cell_setting                     ?',
            '_symmetry.Int_Tables_number                199'
        ]
        exp_cryst1_line_1 = 'CRYST1   77.900   77.900   77.900  90.00  90.00  90.00 I 21 3       24          \n'
        act_cryst1_line_1, exit_1 = make_cryst1_line_from_mmcif(space_group_info_1, False)
        self.assertEqual(exp_cryst1_line_1, act_cryst1_line_1)
        self.assertFalse(exit_1)

        space_group_info_2 = copy.deepcopy(space_group_info_1)
        space_group_info_2[1] = '_cell.length_a       123477.900'
        exp_cryst1_line_2 = ''
        act_cryst1_line_2, exit_2 = make_cryst1_line_from_mmcif(space_group_info_2, False)
        self.assertEqual(exp_cryst1_line_2, act_cryst1_line_2)
        self.assertTrue(exit_2)

        space_group_info_3 = copy.deepcopy(space_group_info_1)
        space_group_info_3[2] = '_cell.length_b                A'
        exp_cryst1_line_3 = ''
        act_cryst1_line_3, exit_3 = make_cryst1_line_from_mmcif(space_group_info_3, False)
        self.assertEqual(exp_cryst1_line_3, act_cryst1_line_3)
        self.assertTrue(exit_3)

        space_group_info_4 = copy.deepcopy(space_group_info_1)
        space_group_info_4[3] = '_cell.length_c 18.5'
        exp_cryst1_line_4 = 'CRYST1   77.900   77.900   18.500  90.00  90.00  90.00 I 21 3       24          \n'
        act_cryst1_line_4, exit_4 = make_cryst1_line_from_mmcif(space_group_info_4, False)
        self.assertEqual(exp_cryst1_line_4, act_cryst1_line_4)
        self.assertFalse(exit_4)

        space_group_info_5 = copy.deepcopy(space_group_info_1)
        space_group_info_5[4] = '_cell.angle_alpha        15.00'
        exp_cryst1_line_5 = 'CRYST1   77.900   77.900   77.900  15.00  90.00  90.00 I 21 3       24          \n'
        act_cryst1_line_5, exit_5 = make_cryst1_line_from_mmcif(space_group_info_5, False)
        self.assertEqual(exp_cryst1_line_5, act_cryst1_line_5)
        self.assertFalse(exit_5)

        space_group_info_6 = copy.deepcopy(space_group_info_1)
        space_group_info_6[5] = '_cell.angle_beta         15.5678'
        exp_cryst1_line_6 = 'CRYST1   77.900   77.900   77.900  90.00  15.57  90.00 I 21 3       24          \n'
        act_cryst1_line_6, exit_6 = make_cryst1_line_from_mmcif(space_group_info_6, False)
        self.assertEqual(exp_cryst1_line_6, act_cryst1_line_6)
        self.assertFalse(exit_6)

        space_group_info_7 = copy.deepcopy(space_group_info_1)
        space_group_info_7[5] = '_cell.angle_gamma         123456'
        exp_cryst1_line_7 = ''
        act_cryst1_line_7, exit_7 = make_cryst1_line_from_mmcif(space_group_info_7, False)
        self.assertEqual(exp_cryst1_line_7, act_cryst1_line_7)
        self.assertTrue(exit_7)

        space_group_info_8 = copy.deepcopy(space_group_info_1)
        space_group_info_8[7] = '_symmetry.space_group_name_H-M \'123456789101\''
        exp_cryst1_line_8 = ''
        act_cryst1_line_8, exit_8 = make_cryst1_line_from_mmcif(space_group_info_8, False)
        self.assertEqual(exp_cryst1_line_8, act_cryst1_line_8)
        self.assertTrue(exit_8)

        space_group_info_9 = copy.deepcopy(space_group_info_1)
        space_group_info_9[8] = '_cell.Z_PDB 12345'
        exp_cryst1_line_9 = ''
        act_cryst1_line_9, exit_9 = make_cryst1_line_from_mmcif(space_group_info_9, False)
        self.assertEqual(exp_cryst1_line_9, act_cryst1_line_9)
        self.assertTrue(exit_9)

        space_group_info_10 = copy.deepcopy(space_group_info_1)
        space_group_info_10[1] = ''
        exp_cryst1_line_10 = ''
        act_cryst1_line_10, exit_10 = make_cryst1_line_from_mmcif(space_group_info_10, False)
        self.assertEqual(exp_cryst1_line_10, act_cryst1_line_10)
        self.assertTrue(exit_10)

    def test_find_disulfides_from_mmcif(self):
        """
        Checks that disulfide bond records in mmCIF format files are parsed as
        expected
        """

        import copy

        # Taken from 2BN1
        disulfide_rec_1 = [
            'loop_',
            '_struct_conn.id',
            '_struct_conn.conn_type_id',
            '_struct_conn.pdbx_leaving_atom_flag',
            '_struct_conn.pdbx_PDB_id',
            '_struct_conn.ptnr1_label_asym_id',
            '_struct_conn.ptnr1_label_comp_id',
            '_struct_conn.ptnr1_label_seq_id',
            '_struct_conn.ptnr1_label_atom_id',
            '_struct_conn.pdbx_ptnr1_label_alt_id',
            '_struct_conn.pdbx_ptnr1_PDB_ins_code',
            '_struct_conn.pdbx_ptnr1_standard_comp_id',
            '_struct_conn.ptnr1_symmetry',
            '_struct_conn.ptnr2_label_asym_id',
            '_struct_conn.ptnr2_label_comp_id',
            '_struct_conn.ptnr2_label_seq_id',
            '_struct_conn.ptnr2_label_atom_id',
            '_struct_conn.pdbx_ptnr2_label_alt_id',
            '_struct_conn.pdbx_ptnr2_PDB_ins_code',
            '_struct_conn.ptnr1_auth_asym_id',
            '_struct_conn.ptnr1_auth_comp_id',
            '_struct_conn.ptnr1_auth_seq_id',
            '_struct_conn.ptnr2_auth_asym_id',
            '_struct_conn.ptnr2_auth_comp_id',
            '_struct_conn.ptnr2_auth_seq_id',
            '_struct_conn.ptnr2_symmetry',
            '_struct_conn.pdbx_ptnr3_label_atom_id',
            '_struct_conn.pdbx_ptnr3_label_seq_id',
            '_struct_conn.pdbx_ptnr3_label_comp_id',
            '_struct_conn.pdbx_ptnr3_label_asym_id',
            '_struct_conn.pdbx_ptnr3_label_alt_id',
            '_struct_conn.pdbx_ptnr3_PDB_ins_code',
            '_struct_conn.details',
            '_struct_conn.pdbx_dist_value',
            '_struct_conn.pdbx_value_order',
            'disulf1 disulf ? ? A CYS 6  SG ? ? ? 1_555 A CYS 11 SG ? ? A CYS 6  A CYS 11 1_555 ? ? ? ? ? ? ? 2.159 ?',
            'disulf2 disulf ? ? A CYS 7  SG ? ? ? 1_555 B CYS 7  SG ? ? A CYS 7  B CYS 7  1_555 ? ? ? ? ? ? ? 2.131 ?',
            'disulf3 disulf ? ? A CYS 20 SG ? ? ? 1_555 B CYS 19 SG ? ? A CYS 20 B CYS 19 1_555 ? ? ? ? ? ? ? 2.018 ?']
        exp_disulfides_1 = {1: [['A', 6, '?'], ['A', 11, '?']],
                            2: [['A', 7, '?'], ['B', 7, '?']],
                            3: [['A', 20, '?'], ['B', 19, '?']]}
        act_disulfides_1, exit_1 = find_disulfides_from_mmcif(disulfide_rec_1, False)
        self.assertDictEqual(act_disulfides_1, exp_disulfides_1)
        self.assertFalse(exit_1)

        act_disulfides_2, exit_2 = find_disulfides_from_mmcif(disulfide_rec_1, True)
        self.assertDictEqual(act_disulfides_2, exp_disulfides_1)
        self.assertTrue(exit_2)


        disulfide_rec_3 = copy.deepcopy(disulfide_rec_1)
        disulfide_rec_3[5] = ''
        act_disulfides_3, exit_3 = find_disulfides_from_mmcif(disulfide_rec_3, False)
        self.assertDictEqual(act_disulfides_3, {})
        self.assertTrue(exit_3)

        disulfide_rec_4 = copy.deepcopy(disulfide_rec_1)
        disulfide_rec_4[-1] = 'disulf3 disulf ? ? A CYS 20 SG ? ? ? 1_555 B CYS'
        act_disulfides_4, exit_4 = find_disulfides_from_mmcif(disulfide_rec_4, False)
        self.assertDictEqual(act_disulfides_4, {})
        self.assertTrue(exit_4)

    def test_parse_seqres_from_mmcif(self):
        """
        Checks that SEQRES records from mmCIF format files are parsed as
        expected
        """

        import copy

        # Taken from 2BN1
        seqres_rec_1 = ['loop_',
                        '_entity_poly_seq.entity_id ',
                        '_entity_poly_seq.num ',
                        '_entity_poly_seq.mon_id ',
                        '_entity_poly_seq.hetero ',
                        '1 1  GLY n ',
                        '1 2  ILE n ',
                        '1 3  VAL n ',
                        '1 4  GLU n ',
                        '1 5  GLN n ',
                        '1 6  CYS n ',
                        '1 7  CYS n ',
                        '1 8  ALA n ',
                        '1 9  SER n ',
                        '1 10 VAL n ',
                        '1 11 CYS n ',
                        '1 12 SER n ',
                        '1 13 LEU n ',
                        '1 14 TYR n ',
                        '1 15 GLN n ',
                        '1 16 LEU n ',
                        '1 17 GLU n ',
                        '1 18 ASN n ',
                        '1 19 TYR n ',
                        '1 20 CYS n ',
                        '1 21 ASN n ',
                        '2 1  PHE n ',
                        '2 2  VAL n ',
                        '2 3  ASN n ',
                        '2 4  GLN n ',
                        '2 5  HIS n ',
                        '2 6  LEU n ',
                        '2 7  CYS n ',
                        '2 8  GLY n ',
                        '2 9  SER n ',
                        '2 10 HIS n ',
                        '2 11 LEU n ',
                        '2 12 VAL n ',
                        '2 13 GLU n ',
                        '2 14 ALA n ',
                        '2 15 LEU n ',
                        '2 16 TYR n ',
                        '2 17 LEU n ',
                        '2 18 VAL n ',
                        '2 19 CYS n ',
                        '2 20 GLY n ',
                        '2 21 GLU n ',
                        '2 22 ARG n ',
                        '2 23 GLY n ',
                        '2 24 PHE n ',
                        '2 25 PHE n ',
                        '2 26 TYR n ',
                        '2 27 THR n ',
                        '2 28 PRO n ',
                        '2 29 LYS n ',
                        '2 30 ALA n ']
        exp_seq_1 = [
            'GLY', 'ILE', 'VAL', 'GLU', 'GLN', 'CYS', 'ALA', 'SER', 'LEU',
            'TYR', 'ASN', 'PHE', 'HIS', 'ARG', 'THR', 'PRO', 'LYS'
        ]
        act_seq_1, exit_1 = parse_seqres_from_mmcif(seqres_rec_1, False)
        self.assertEqual(exp_seq_1, act_seq_1)
        self.assertFalse(exit_1)

        act_seq_2, exit_2 = parse_seqres_from_mmcif(seqres_rec_1, True)
        self.assertEqual(exp_seq_1, act_seq_2)
        self.assertTrue(exit_2)

        seqres_rec_3 = copy.deepcopy(seqres_rec_1)
        seqres_rec_3[3] = ''
        exp_seq_3 = []
        act_seq_3, exit_3 = parse_seqres_from_mmcif(seqres_rec_3, False)
        self.assertEqual(exp_seq_3, act_seq_3)
        self.assertTrue(exit_3)

    def test_parse_atomrec_from_mmcif(self):
        """
        Checks that ATOM/HETATM records in mmCIF format files are parsed as
        expected
        """

        import copy

        # Lines taken from 2BN1 and 2QTZ
        atom_lines_1 = [
            'loop_',
            '_atom_site.group_PDB ',
            '_atom_site.id ',
            '_atom_site.type_symbol ',
            '_atom_site.label_atom_id ',
            '_atom_site.label_alt_id ',
            '_atom_site.label_comp_id ',
            '_atom_site.label_asym_id ',
            '_atom_site.label_entity_id ',
            '_atom_site.label_seq_id ',
            '_atom_site.pdbx_PDB_ins_code ',
            '_atom_site.Cartn_x ',
            '_atom_site.Cartn_y ',
            '_atom_site.Cartn_z ',
            '_atom_site.occupancy ',
            '_atom_site.B_iso_or_equiv ',
            '_atom_site.pdbx_formal_charge ',
            '_atom_site.auth_seq_id ',
            '_atom_site.auth_comp_id ',
            '_atom_site.auth_asym_id ',
            '_atom_site.auth_atom_id ',
            '_atom_site.pdbx_PDB_model_num ',
            'ATOM   2    C CA . GLY A 1 1  ? 13.602 45.768 30.728 1.00 19.08 ? 1    GLY A CA  1 ',
            'HETATM 438  O O  . HOH C 3 .  ? 13.128 46.298 34.342 1.00 33.10 ? 2001 HOH A O   1 ',
            'HETATM 3546 P PA . FAD B 2 .  ? 29.643 29.700 51.402 1.00 20.52 ? 700  FAD A PA  1 ']

        exp_atom_rec_1 = []

        exp_atom_1 = atom()
        exp_atom_1.lineID = 'ATOM'
        exp_atom_1.atomNum = 2
        exp_atom_1.atomType = 'CA'
        exp_atom_1.conformer = ''
        exp_atom_1.resiType = 'GLY'
        exp_atom_1.chainID = 'A'
        exp_atom_1.resiNum = 1
        exp_atom_1.insCode = ''
        exp_atom_1.xyzCoords = [[13.602], [45.768], [30.728]]
        exp_atom_1.occupancy = 1.0
        exp_atom_1.bFactor = 19.08
        exp_atom_1.element = 'C'
        exp_atom_1.charge = ''
        exp_atom_1.origResiNum = '1'
        exp_atom_1.origResiType = 'GLY'
        exp_atom_1.origChainID = 'A'
        exp_atom_1.origAtomType = 'CA'
        exp_atom_1.pd = None
        exp_atom_1.avrg_bf = None
        exp_atom_1.bd = None
        exp_atom_rec_1.append(exp_atom_1)

        exp_atom_2 = atom()
        exp_atom_2.lineID = 'HETATM'
        exp_atom_2.atomNum = 438
        exp_atom_2.atomType = 'O'
        exp_atom_2.conformer = ''
        exp_atom_2.resiType = 'HOH'
        exp_atom_2.chainID = 'A'
        exp_atom_2.resiNum = 2001
        exp_atom_2.insCode = ''
        exp_atom_2.xyzCoords = [[13.128], [46.298], [34.342]]
        exp_atom_2.occupancy = 1.0
        exp_atom_2.bFactor = 33.1
        exp_atom_2.element = 'O'
        exp_atom_2.charge = ''
        exp_atom_2.origResiNum = ''
        exp_atom_2.origResiType = 'HOH'
        exp_atom_2.origChainID = 'C'
        exp_atom_2.origAtomType = 'O'
        exp_atom_2.pd = None
        exp_atom_2.avrg_bf = None
        exp_atom_2.bd = None
        exp_atom_rec_1.append(exp_atom_2)

        exp_atom_3 = atom()
        exp_atom_3.lineID = 'HETATM'
        exp_atom_3.atomNum = 3546
        exp_atom_3.atomType = 'PA'
        exp_atom_3.conformer = ''
        exp_atom_3.resiType = 'FAD'
        exp_atom_3.chainID = 'A'
        exp_atom_3.resiNum = 700
        exp_atom_3.insCode = ''
        exp_atom_3.xyzCoords = [[29.643], [29.7], [51.402]]
        exp_atom_3.occupancy = 1.0
        exp_atom_3.bFactor = 20.52
        exp_atom_3.element = 'P'
        exp_atom_3.charge = ''
        exp_atom_3.origResiNum = ''
        exp_atom_3.origResiType = 'FAD'
        exp_atom_3.origChainID = 'B'
        exp_atom_3.origAtomType = 'PA'
        exp_atom_3.pd = None
        exp_atom_3.avrg_bf = None
        exp_atom_3.bd = None
        exp_atom_rec_1.append(exp_atom_3)

        act_atom_rec_1, exit_1 = parse_atom_rec_from_mmcif(atom_lines_1, False)
        for index in range(len(act_atom_rec_1)):
            self.assertEqual(exp_atom_rec_1[index], act_atom_rec_1[index])
        self.assertFalse(exit_1)

        act_atom_rec_2, exit_2 = parse_atom_rec_from_mmcif(atom_lines_1, True)
        for index in range(len(act_atom_rec_2)):
            self.assertEqual(exp_atom_rec_1[index], act_atom_rec_2[index])
        self.assertTrue(exit_2)

        atom_lines_3 = copy.deepcopy(atom_lines_1)
        atom_lines_3[-1] = 'HETATM X P PA . FAD B 2 .  ? 29.643 29.700 51.402 1.00 20.52 ? 700  FAD A PA  1 '
        exp_atom_rec_3 = []
        act_atom_rec_3, exit_3 = parse_atom_rec_from_mmcif(atom_lines_3, False)
        self.assertEqual(exp_atom_rec_3, act_atom_rec_3)
        self.assertTrue(exit_3)

        atom_lines_4 = copy.deepcopy(atom_lines_1)
        atom_lines_4[-1] = 'HETATM 3546 P PA . FAD B 2 .  ? 29.643 29.700 51.402 1.00 20.52 ? X  FAD A PA  1 '
        exp_atom_rec_4 = []
        act_atom_rec_4, exit_4 = parse_atom_rec_from_mmcif(atom_lines_4, False)
        self.assertEqual(exp_atom_rec_4, act_atom_rec_4)
        self.assertTrue(exit_4)

        atom_lines_5 = copy.deepcopy(atom_lines_1)
        atom_lines_5[-1] = 'HETATM 3546 P PA . FAD B 2 .  ? 29.643 X 51.402 1.00 20.52 ? 700  FAD A PA  1 '
        exp_atom_rec_5 = []
        act_atom_rec_5, exit_5 = parse_atom_rec_from_mmcif(atom_lines_5, False)
        self.assertEqual(exp_atom_rec_5, act_atom_rec_5)
        self.assertTrue(exit_5)

        atom_lines_6 = copy.deepcopy(atom_lines_1)
        atom_lines_6[-1] = 'HETATM 3546 P PA . FAD B 2 .  ? 29.643 29.700 51.402 X 20.52 ? 700  FAD A PA  1 '
        exp_atom_rec_6 = []
        act_atom_rec_6, exit_6 = parse_atom_rec_from_mmcif(atom_lines_6, False)
        self.assertEqual(exp_atom_rec_6, act_atom_rec_6)
        self.assertTrue(exit_6)

        atom_lines_7 = copy.deepcopy(atom_lines_1)
        atom_lines_7[-1] = 'HETATM 3546 P PA . FAD B 2 .  ? 29.643 29.700 51.402 1.00 X ? 700  FAD A PA  1 '
        exp_atom_rec_7 = []
        act_atom_rec_7, exit_7 = parse_atom_rec_from_mmcif(atom_lines_7, False)
        self.assertEqual(exp_atom_rec_7, act_atom_rec_7)
        self.assertTrue(exit_7)

        atom_lines_8 = copy.deepcopy(atom_lines_1)
        atom_lines_8[2] = ''
        exp_atom_rec_8 = []
        act_atom_rec_8, exit_8 = parse_atom_rec_from_mmcif(atom_lines_8, False)
        self.assertEqual(exp_atom_rec_8, act_atom_rec_8)
        self.assertTrue(exit_8)

    def test_write_mmcif_file(self):
        """
        Tests that atom objects are correctly converted into mmCIF format
        """

        import copy
        import pandas as pd


        atoms_list = []
        # Taken from 2BN1 and 2QTZ
        new_atom_1 = atom()
        new_atom_1.lineID = 'ATOM'
        new_atom_1.atomNum = 2
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
        new_atom_1.origResiNum = '1'
        new_atom_1.origResiType = 'GLY'
        new_atom_1.origChainID = 'A'
        new_atom_1.origAtomType = 'CA'
        new_atom_1.pd = 19
        new_atom_1.avrg_bf = 26.2
        new_atom_1.bd = 1.3
        atoms_list.append(new_atom_1)

        new_atom_2 = atom()
        new_atom_2.lineID = 'HETATM'
        new_atom_2.atomNum = 438
        new_atom_2.atomType = 'O'
        new_atom_2.conformer = ''
        new_atom_2.resiType = 'HOH'
        new_atom_2.chainID = 'B'
        new_atom_2.resiNum = 2001
        new_atom_2.insCode = ''
        new_atom_2.xyzCoords = [[13.128], [46.298], [34.342]]
        new_atom_2.occupancy = 1.0
        new_atom_2.bFactor = 33.1
        new_atom_2.element = 'O'
        new_atom_2.charge = ''
        new_atom_2.origResiNum = ''
        new_atom_2.origResiType = 'HOH'
        new_atom_2.origChainID = 'C'
        new_atom_2.origAtomType = 'O'
        new_atom_2.pd = 30
        new_atom_2.avrg_bf = 25.6
        new_atom_2.bd = 0.97
        atoms_list.append(new_atom_2)

        new_atom_3 = atom()
        new_atom_3.lineID = 'HETATM'
        new_atom_3.atomNum = 3546
        new_atom_3.atomType = 'PA'
        new_atom_3.conformer = 'X'
        new_atom_3.resiType = 'FAD'
        new_atom_3.chainID = 'B'
        new_atom_3.resiNum = 700
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
        atoms_list.append(new_atom_3)

        exp_cif_lines = [
            '#',
            'loop_',
            '_atom_site.group_PDB',
            '_atom_site.id',
            '_atom_site.type_symbol',
            '_atom_site.label_atom_id',
            '_atom_site.label_alt_id',
            '_atom_site.label_comp_id',
            '_atom_site.label_asym_id',
            '_atom_site.label_seq_id',
            '_atom_site.pdbx_PDB_ins_code',
            '_atom_site.Cartn_x',
            '_atom_site.Cartn_y',
            '_atom_site.Cartn_z',
            '_atom_site.occupancy',
            '_atom_site.B_iso_or_equiv',
            '_atom_site.B_Damage',
            '_atom_site.packing_density',
            '_atom_site.average_B_factor',
            '_atom_site.pdbx_formal_charge',
            '_atom_site.auth_seq_id',
            '_atom_site.auth_comp_id',
            '_atom_site.auth_asym_id',
            '_atom_site.auth_atom_id',
            'ATOM   2    C CA . GLY A 1 ? 13.602 45.768 30.728 1.0 19.08 1.3  19 26.2 ? 1    GLY A CA',
            'HETATM 438  O O  . HOH C . ? 13.128 46.298 34.342 1.0 33.1  0.97 30 25.6 ? 2001 HOH B O ',
            'HETATM 3546 P PA X FAD B . ? 29.643 29.7   51.402 1.0 20.52 2.0  26 20.5 ? 700  FAD B PA',
            '#',
            ''
        ]

        empty_df = pd.DataFrame({})
        output = generate_output_files('tests/temp_files/test', '', empty_df)
        output.write_output_cif(atoms_list)
        with open('tests/temp_files/test_BDamage.cif', 'r') as f:
            act_cif_lines = f.read().split('\n')

        self.assertEqual(exp_cif_lines, act_cif_lines)
