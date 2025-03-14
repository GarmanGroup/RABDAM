
# RABDAM
# Copyright (C) 2025 Garman Group, University of Oxford

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


# python -m unittest tests/test_model_filtering.py

import unittest
from rabdam.Subroutines.CalculateBDamage import run_rabdam
from rabdam.Subroutines.PDBCUR import check_for_protein


class TestClass(unittest.TestCase):

    def test_check_for_protein(self):
        """
        Checks that recognition of protein chains works
        """

        import os
        import requests
        import shutil

        if not os.path.isdir('tests/temp_files/'):
            os.mkdir('tests/temp_files/')

        # Check DNA structure
        cif_text = requests.get('https://files.rcsb.org/view/6QT2.cif')
        with open('tests/temp_files/6QT2.cif', 'w') as f:
            f.write(cif_text.text)
        contains_protein_6qt2 = check_for_protein('tests/temp_files/6QT2.cif')
        self.assertFalse(contains_protein_6qt2)

        # Check protein structure
        cif_text = requests.get('https://files.rcsb.org/view/2BN1.cif')
        with open('tests/temp_files/2BN1.cif', 'w') as f:
            f.write(cif_text.text)
        contains_protein_2bn1 = check_for_protein('tests/temp_files/2BN1.cif')
        self.assertTrue(contains_protein_2bn1)

        # Check protein/NA complex structure
        cif_text = requests.get('https://files.rcsb.org/view/5EEU.cif')
        with open('tests/temp_files/5EEU.cif', 'w') as f:
            f.write(cif_text.text)
        contains_protein_5eeu = check_for_protein('tests/temp_files/5EEU.cif')
        self.assertTrue(contains_protein_5eeu)

        shutil.rmtree('tests/temp_files/')

    def test_model_filtering(self):
        """
        Checks that models that don't meet requirements for Bnet calculation set
        out in Shelley & Garman, 2022 are successfully filtered out
        """

        import os
        import requests
        import shutil

        pdb_codes = [
            '1AWQ',  # Rfree > 0.4
            '1H3Y',  # Resolution > 3.5
            '101M',  # Temperature > 120K
            '2BLZ',  # Sub-1 occupancy asp/glu
            '6QT2',  # DNA
            '1A38',  # Per-residue B-factors
            '6Q8T',  # < 20 Asp + Glu O atoms
        ]

        if not os.path.isdir('tests/temp_files/'):
            os.mkdir('tests/temp_files/')

        cwd = os.getcwd()
        for code in pdb_codes:
            # Downloads cif file
            cif_text = requests.get('https://files.rcsb.org/view/%s.cif' % code)
            with open('tests/temp_files/%s.cif' % code, 'w') as f:
                f.write(cif_text.text)

            # Run RABDAM with filter
            rabdam_filter = run_rabdam(
                pathToInput='%s/tests/temp_files/%s.cif' % (cwd, code),
                outputDir='%s/tests/temp_files/' % cwd,
                batchRun=True,
                overwrite=True,
                outFiles="all",
                filterInput=True,
                temperature=None,
                resolution=None,
                PDT=7,
                windowSize=0.02,
                HETATM=False,
                removeAtoms=[],
                addAtoms=[],
                highlightAtoms=[],
                createOrigpdb=False,
                createAUpdb=False,
                createUCpdb=False,
                createAUCpdb=False,
                createTApdb=False
                )
            success = rabdam_filter.rabdam_dataframe(test=True)
            self.assertFalse(success)

            # Run RABDAM without filter
            rabdam_no_filter = run_rabdam(
                pathToInput='%s/tests/temp_files/%s.cif' % (cwd, code),
                outputDir='%s/tests/temp_files/' % cwd,
                batchRun=True,
                overwrite=True,
                outFiles="all",
                filterInput=False,
                temperature=None,
                resolution=None,
                PDT=7,
                windowSize=0.02,
                HETATM=False,
                removeAtoms=[],
                addAtoms=[],
                highlightAtoms=[],
                createOrigpdb=False,
                createAUpdb=False,
                createUCpdb=False,
                createAUCpdb=False,
                createTApdb=False
                )
            success = rabdam_no_filter.rabdam_dataframe(test=True)
            if code == '6QT2':  # Will fail for DNA because all non-protein
            # atoms will be filtered out
                self.assertFalse(success)
            else:
                self.assertTrue(success)

            # Check resolution and temperature overrides in input file
            if code == '1H3Y':
                # Downloads cif file
                cif_text = requests.get('https://files.rcsb.org/view/%s.cif' % code)
                with open('tests/temp_files/%s.cif' % code, 'w') as f:
                    f.write(cif_text.text)

                # Runs RABDAM
                rabdam = run_rabdam(
                    pathToInput='%s/tests/temp_files/%s.cif' % (os.getcwd(), code),
                    outputDir='%s/tests/temp_files/' % os.getcwd(),
                    batchRun=True,
                    overwrite=True,
                    outFiles="all",
                    filterInput=True,
                    temperature=None,
                    resolution=1.5,
                    PDT=7,
                    windowSize=0.02,
                    HETATM=False,
                    removeAtoms=[],
                    addAtoms=[],
                    highlightAtoms=[],
                    createOrigpdb=False,
                    createAUpdb=False,
                    createUCpdb=False,
                    createAUCpdb=False,
                    createTApdb=False
                )
                success = rabdam.rabdam_dataframe(test=True)
                self.assertTrue(success)

            # Check temperature override in input file
            elif code == '101M':
                # Downloads cif file
                cif_text = requests.get('https://files.rcsb.org/view/%s.cif' % code)
                with open('tests/temp_files/%s.cif' % code, 'w') as f:
                    f.write(cif_text.text)

                # Runs RABDAM
                rabdam = run_rabdam(
                    pathToInput='%s/tests/temp_files/%s.cif' % (os.getcwd(), code),
                    outputDir='%s/tests/temp_files/' % os.getcwd(),
                    batchRun=True,
                    overwrite=True,
                    outFiles="all",
                    filterInput=True,
                    temperature=100,
                    resolution=None,
                    PDT=7,
                    windowSize=0.02,
                    HETATM=False,
                    removeAtoms=[],
                    addAtoms=[],
                    highlightAtoms=[],
                    createOrigpdb=False,
                    createAUpdb=False,
                    createUCpdb=False,
                    createAUCpdb=False,
                    createTApdb=False
                )
                success = rabdam.rabdam_dataframe(test=True)
                self.assertTrue(success)

        shutil.rmtree('tests/temp_files/')
        
            
            

        