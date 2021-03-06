
# RABDAM
# Copyright (C) 2020 Garman Group, University of Oxford

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

# An outer layer to the pipeline scripts. Depending upon the flags specified
# in the command line input, this script will run either the complete / a
# subsection of the pipeline.

# python -m unittest tests/test_parse_input_file.py

import unittest
from rabdam.rabdam import (
    ArgumentError, FileDoesNotExistError, parse_command_line_arguments,
    parse_input_file_arguments
)

class TestClass(unittest.TestCase):

    def test_parse_command_line_arguments(self):
        """
        Checks that command line arguments are parsed as expected
        """

        import argparse

        command_line_1 = []
        self.assertRaises(
            SystemExit, parse_command_line_arguments, command_line_1, True
        )

        command_line_3 = ['-i', 'tests/test_files/example input.txt']
        act_parsed_3 = parse_command_line_arguments(command_line_3, True)
        exp_parsed_3 = argparse.Namespace(
            dependencies=False, input='tests/test_files/example input.txt',
            pdb_or_mmcif_file=None, output=['csv', 'pdb', 'cif', 'kde', 'bnet',
            'summary'], run=None
        )
        self.assertEqual(act_parsed_3, exp_parsed_3)

        command_line_4 = ['-i', 'does_not_exist.txt']
        self.assertRaises(
            FileDoesNotExistError, parse_command_line_arguments, command_line_4,
            True
        )

        command_line_5 = ['-f', '2bn1']
        act_parsed_5 = parse_command_line_arguments(command_line_5, True)
        exp_parsed_5 = argparse.Namespace(
            dependencies=False, input=None, pdb_or_mmcif_file=['2bn1'],
            output=['csv', 'pdb', 'cif', 'kde', 'bnet', 'summary'], run=None
        )
        self.assertEqual(act_parsed_5, exp_parsed_5)

        command_line_6 = ['-f', '5eeu', '5eev', '5eew']
        act_parsed_6 = parse_command_line_arguments(command_line_6, True)
        exp_parsed_6 = argparse.Namespace(
            dependencies=False, input=None, pdb_or_mmcif_file=['5eeu', '5eev',
            '5eew'], output=['csv', 'pdb', 'cif', 'kde', 'bnet', 'summary'],
            run=None
        )
        self.assertEqual(act_parsed_6, exp_parsed_6)

        command_line_7 = ['-o', 'csv', 'cif', 'pdb']
        self.assertRaises(
            SystemExit, parse_command_line_arguments, command_line_7, True
        )

        command_line_8 = ['-o', 'csv', 'cif', 'pdb', '-f', '2bn1']
        act_parsed_8 = parse_command_line_arguments(command_line_8, True)
        exp_parsed_8 = argparse.Namespace(
            dependencies=False, input=None, pdb_or_mmcif_file=['2bn1'],
            output=['csv', 'cif', 'pdb'], run=None
        )
        self.assertEqual(act_parsed_8, exp_parsed_8)

        command_line_9 = ['-o', 'csv', 'CIF', '-f', '2bn1']
        act_parsed_9 = parse_command_line_arguments(command_line_9, True)
        exp_parsed_9 = argparse.Namespace(
            dependencies=False, input=None, pdb_or_mmcif_file=['2bn1'],
            output=['csv', 'cif'], run=None
        )
        self.assertEqual(act_parsed_9, exp_parsed_9)

        command_line_10 = ['-o', 'csv', 'cfi', '-f', '2bn1']
        self.assertRaises(
            ArgumentError, parse_command_line_arguments, command_line_10, True
        )

        command_line_11 = ['-r', 'df', '-f', '2bn1']
        act_parsed_11 = parse_command_line_arguments(command_line_11, True)
        exp_parsed_11 = argparse.Namespace(
            dependencies=False, input=None, pdb_or_mmcif_file=['2bn1'],
            output=['csv', 'pdb', 'cif', 'kde', 'bnet', 'summary'], run='df'
        )
        self.assertEqual(act_parsed_11, exp_parsed_11)

        command_line_12 = ['-r', 'anaylsis', '-f', '2bn1']
        self.assertRaises(
            ArgumentError, parse_command_line_arguments, command_line_12, True
        )

    def test_parse_input_file(self):
        """
        Checks that program arguments listed in input file are parsed as expected
        """

        import os

        cwd = os.getcwd()

        input_1 = []
        act_output_1 = parse_input_file_arguments(input_1)
        exp_output_1 = {'outputDir': cwd,
                        'batchRun': False,
                        'overwrite': False,
                        'PDT': 7.0,
                        'windowSize': 0.02,
                        'protOrNA': 'protein',
                        'HETATM': False,
                        'addAtoms': [],
                        'removeAtoms': [],
                        'highlightAtoms': [],
                        'createOrigpdb': False,
                        'createAUpdb': False,
                        'createUCpdb': False,
                        'createAUCpdb': False,
                        'createTApdb': False}
        self.assertDictEqual(act_output_1, exp_output_1)

        input_2 = [
            'outputdir=tests/', 'batchcontinue=y', 'overwrite=n', 'pdt=4',
            'windowsize=0.5', 'proteinornucleicacid=protein', 'hetatm=keep',
            'removeatoms=2;5;7-9;HOH', 'addatoms=NA;FE;45-49;3;FAD',
            'highlightatoms=1-6', 'createorigpdb=YES', 'createaupdb=No',
            'createucpdb=true', 'createaucpdb=False', 'createtapdb=True']
        act_output_2 = parse_input_file_arguments(input_2)
        exp_output_2 = {'outputDir': 'tests/',
                        'batchRun': True,
                        'overwrite': False,
                        'PDT': 4.0,
                        'windowSize': 0.5,
                        'protOrNA': 'protein',
                        'HETATM': True,
                        'removeAtoms': ['2', '5', '7', '8', '9', 'HOH'],
                        'addAtoms': ['NA', 'FE', '45', '46', '47', '48', '49', '3', 'FAD'],
                        'highlightAtoms': ['1', '2', '3', '4', '5', '6'],
                        'createOrigpdb': True,
                        'createAUpdb': False,
                        'createUCpdb': True,
                        'createAUCpdb': False,
                        'createTApdb': True}
        self.assertDictEqual(act_output_2, exp_output_2)

        input_3 = ['outputdir=does_not_exist/']
        self.assertRaises(FileDoesNotExistError, parse_input_file_arguments, input_3)

        input_4 = ['batchcontinue=xyz']
        self.assertRaises(ArgumentError, parse_input_file_arguments, input_4)

        input_5 = ['overwrite=xyz']
        self.assertRaises(ArgumentError, parse_input_file_arguments, input_5)

        input_6 = ['pdt=xyz']
        self.assertRaises(ArgumentError, parse_input_file_arguments, input_6)

        input_7 = ['pdt=1.5']
        self.assertRaises(ArgumentError, parse_input_file_arguments, input_7)

        input_8 = ['windowSize=xyz']
        self.assertRaises(ArgumentError, parse_input_file_arguments, input_8)

        input_9 = ['windowSize=1.5']
        self.assertRaises(ArgumentError, parse_input_file_arguments, input_9)

        input_10 = ['hetatm=True']
        self.assertRaises(ArgumentError, parse_input_file_arguments, input_10)

        input_11 = ['removeAtoms=A-C']
        self.assertRaises(ArgumentError, parse_input_file_arguments, input_11)

        input_12 = ['addAtoms=1-2-3']
        self.assertRaises(ArgumentError, parse_input_file_arguments, input_12)

        input_13 = ['highlightAtoms=1-A']
        self.assertRaises(ArgumentError, parse_input_file_arguments, input_13)

        input_14 = ['removeAtoms=A-C']
        self.assertRaises(ArgumentError, parse_input_file_arguments, input_14)

        input_15 = ['createOrigpdb=Flase']
        self.assertRaises(ArgumentError, parse_input_file_arguments, input_15)

        input_16 = ['createAUpdb=Yse']
        self.assertRaises(ArgumentError, parse_input_file_arguments, input_16)

        input_17 = ['createUCpdb=1-6']
        self.assertRaises(ArgumentError, parse_input_file_arguments, input_17)

        input_18 = ['createAUCpdb=blah']
        self.assertRaises(ArgumentError, parse_input_file_arguments, input_18)

        input_19 = ['createTApdb=np.inf']
        self.assertRaises(ArgumentError, parse_input_file_arguments, input_19)

        input_20 = ['unrecognised=1']
        self.assertRaises(ArgumentError, parse_input_file_arguments, input_20)
