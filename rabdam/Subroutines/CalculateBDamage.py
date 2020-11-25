
# RABDAM
# Copyright (C) 2018 Garman Group, University of Oxford

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


class rabdam(object):
    def __init__(self, pathToInput, outputDir, batchRun, overwrite, PDT,
                 windowSize, protOrNA, HETATM, addAtoms, removeAtoms,
                 highlightAtoms, createOrigpdb, createAUpdb, createUCpdb,
                 createAUCpdb, createTApdb):
        self.pathToInput = pathToInput
        self.outputDir = outputDir
        self.batchRun = batchRun
        self.overwrite = overwrite
        self.PDT = PDT
        self.windowSize = windowSize
        self.protOrNA = protOrNA
        self.HETATM = HETATM
        self.addAtoms = addAtoms
        self.removeAtoms = removeAtoms
        self.highlightAtoms = highlightAtoms
        self.createOrigpdb = createOrigpdb
        self.createAUpdb = createAUpdb
        self.createUCpdb = createUCpdb
        self.createAUCpdb = createAUCpdb
        self.createTApdb = createTApdb

    def rabdam_dataframe(self, run):
        """
        Calculates BDamage for selected atoms within input PDB file and
        writes output to DataFrame.
        """

        prompt = '> '
        import sys
        import os
        import shutil
        import requests
        import numpy as np
        import pickle

        if sys.version_info[0] < 3:
            user_input = raw_input
        else:
            user_input = input

        if __name__ == 'Subroutines.CalculateBDamage':
            from Subroutines.PDBCUR import (
                parse_mmcif_file, parse_pdb_file, clean_atom_rec, genPDBCURinputs,
                runPDBCUR
                )
            from Subroutines.parsePDB import (
                download_mmcif, copy_input, full_atom_list, b_damage_atom_list
                )
            from Subroutines.translateUnitCell import (
                convertToCartesian, translateUnitCell, extract_unit_cell_params
                )
            from Subroutines.trimUnitCellAssembly import (
                getAUparams, convertParams, trimAtoms
                )
            from Subroutines.makeDataFrame import (
                convert_array_to_atom_list, makePDB, writeDataFrame
            )
            from Subroutines.BDamage import (
                get_xyz_from_objects, calc_packing_density,
                write_pckg_dens_to_atoms, calcBDam
                )
        else:
            from rabdam.Subroutines.PDBCUR import (
                parse_mmcif_file, parse_pdb_file, clean_atom_rec, genPDBCURinputs,
                runPDBCUR
                )
            from rabdam.Subroutines.parsePDB import (
                download_mmcif, copy_input, full_atom_list, b_damage_atom_list
                )
            from rabdam.Subroutines.translateUnitCell import (
                convertToCartesian, translateUnitCell, extract_unit_cell_params
                )
            from rabdam.Subroutines.trimUnitCellAssembly import (
                getAUparams, convertParams, trimAtoms
                )
            from rabdam.Subroutines.makeDataFrame import (
                convert_array_to_atom_list, makePDB, writeDataFrame
            )
            from rabdam.Subroutines.BDamage import (
                get_xyz_from_objects, calc_packing_density,
                write_pckg_dens_to_atoms, calcBDam
                )

        if run == 'rabdam':
            print('**************************** RABDAM ****************************\n')
        elif run == 'rabdam_dataframe':
            print('*********************** RABDAM DATAFRAME ***********************\n')

        print('\n****************************************************************\n'
              '***************** Program to calculate BDamage *****************\n'
              '****************************************************************\n')

        print('****************************************************************\n'
              '************************* Input Section ************************\n')
        # Prints the values of the program options to be used in the current
        # RABDAM run
        print('Calculating BDamage for %s' % self.pathToInput)
        print('Writing output files to %s' % self.outputDir)
        if self.PDT == 7:
            print('Using default packing density threshold of 7 Angstroms')
        else:
            print('Packing density threshold defined by user as %s Angstroms' % self.PDT)
        if self.windowSize == 0.02:
            print('Using default window size of 2%')
        else:
            print('Window size defined by user as %s%%' % (self.windowSize*100))
        if self.HETATM is True:
            print('Keeping HETATM')
        elif self.HETATM is False:
            print('Removing HETATM (default)')
        if self.protOrNA == 'protein':
            print('Retaining protein atoms, discarding nucleic acid atoms (default)')
        elif self.protOrNA in ['nucleicacid', 'na']:
            print('Retaining nucleic acid atoms, discarding protein atoms')
        if len(self.removeAtoms) == 0:
            print('No atoms to be removed (default)')
        else:
            remove_atoms_string = ''
            for value in self.removeAtoms:
                remove_atoms_string = remove_atoms_string + value + ', '
            print('Atoms to be removed: %s' % remove_atoms_string[:-2])
        if len(self.addAtoms) == 0:
            print('No atoms to be added (default)')
        else:
            add_atoms_string = ''
            for value in self.addAtoms:
                add_atoms_string = add_atoms_string + value + ', '
            print('Atoms to be added: %s' % add_atoms_string[:-2])

        print('\n********************* End of Input Section *********************\n'
              '****************************************************************\n')

        # Changes directory to the specified location for the output 'Logfiles'
        # directory. The default location is the current working directory.
        cwd = os.getcwd()
        os.chdir(self.outputDir)

        print('****************************************************************\n'
              '********************** Process PDB Section *********************\n')
        # Creates a new directory named after the input PDB file (after
        # checking that this directory does not already exist) in the
        # 'Logfiles' directory. Then saves a copy of the input PDB file to the
        # new directory.

        # If 4 digit PDB accession code has been supplied:
        if len(self.pathToInput) == 4:
            print('Accession code supplied')
            PDBcode = self.pathToInput.upper()
            window_name = str(100*self.windowSize).replace('.', '_')
            pdt_name = str(self.PDT).replace('.', '_')
            PDBdirectory = 'Logfiles/%s_window_%s_pdt_%s/' % (PDBcode, window_name, pdt_name)
            file_name_start = '%s%s' % (PDBdirectory, PDBcode)
            pathToInput = '%s%s.cif' % (PDBdirectory, PDBcode)

            # Checks whether accession code exists - if not, exit program
            # with error message
            mmcif_url = 'https://files.rcsb.org/view/%s.cif' % PDBcode
            header = requests.get(mmcif_url)
            if header.status_code != 200:
                print('\n\nERROR: Failed to download %s mmCIF file with '
                      'accession code %s:\nCheck that a structure with this '
                      'accession code exists.' % (url[-4:], PDBcode))
                if self.batchRun is False:
                    sys.exit()
                elif self.batchRun is True:
                    return

            # If directory with same name as PDBdirectory already exists in
            # 'Logfiles' directory, user input is requested ('Do you want to
            # overwrite the existing folder?'):
            # yes = old PDBdirectory is deleted, new PDBdirectory is created
            #       and copy of the mmCIF file is downloaded from the RSCB PDB
            #       website and saved to the new directory
            # no = old PDBdirectory is retained, exit program
            # If it doesn't already exist, new PDBdirectory is created and
            # copy of the mmCIF file is downloaded from the RSCB PDB website
            # and saved to the new directory.
            if os.path.isdir(PDBdirectory):
                print('\nFolder %s already exists locally at %s' % (
                    PDBcode, PDBdirectory
                ))
                print('Do you want to overwrite the existing folder?\n'
                      '--USER INPUT-- type your choice and press RETURN\n'
                      'yes = overwrite this folder\n'
                      'no = do not overwrite this folder\n')
                owChoice = None
                while owChoice not in ['yes', 'no', 'y', 'n']:
                    if self.overwrite is True:
                        owChoice = 'yes'
                    elif self.overwrite is False and self.batchRun is True:
                        owChoice = 'no'
                    elif self.overwrite is False and self.batchRun is False:
                        owChoice = user_input(prompt).lower()

                    if owChoice == 'yes' or owChoice == 'y':
                        print('\nOverwriting existing folder')
                        shutil.rmtree(PDBdirectory)
                        download_mmcif(PDBcode, PDBdirectory, pathToInput)
                        break
                    elif owChoice == 'no' or owChoice == 'n':
                        print('\nKeeping original folder\n')
                        if self.batchRun is False:
                            sys.exit()
                        elif self.batchRun is True:
                            return
                        break
                    else:
                        print('Unrecognised input - please answer "yes" or "no"')
            else:
                download_mmcif(PDBcode, PDBdirectory, pathToInput)

            # Checks that mmCIF file has been successfully downloaded and saved
            # to the 'Logfiles' directory
            if not os.path.exists(pathToInput):
                shutil.rmtree('%s' % PDBdirectory)
                print('\n\nERROR: Failed to download and save %s mmCIF file - '
                      'cause unknown' % PDBcode)
                if self.batchRun is False:
                    sys.exit()
                elif self.batchRun is True:
                    return

        # If filepath to PDB / mmCIF file has been supplied:
        else:
            # Changes directory to allow input PDB / mmCIF file to be read from
            # any provided file path. If directory with same name as
            # PDBdirectory already exists in 'Logfiles' directory, user input
            # is requested ('Do you want to overwrite the existing folder?'):
            # yes = old PDBdirectory is deleted, new PDBdirectory is created
            #       and copy of input PDB / mmCIF file is saved to the new
            #       directory
            # no = old PDBdirectory is retained, exit program
            # If it doesn't already exist, new PDBdirectory is created and copy
            # of input PDB / mmCIF file is saved to the new directory.
            pathToInput = '/' + self.pathToInput.replace('\\', '/').strip('/')
            splitPath = pathToInput.split('/')
            disk = '%s/' % splitPath[0]
            os.chdir('/')
            os.chdir(disk)
            if not os.path.exists(pathToInput):
                print('\n\nERROR: Supplied filepath not recognised')
                if self.batchRun is False:
                    sys.exit()
                elif self.batchRun is True:
                    return
            os.chdir(self.outputDir)

            if pathToInput[-4:] not in ['.pdb', '.cif']:
                print('\n\nERROR: Supplied input filepath is not a .pdb or '
                      '.cif file')
                if self.batchRun is False:
                    sys.exit()
                elif self.batchRun is True:
                    return

            print('Filepath to %s file supplied\n' % pathToInput[-4:])
            splitFilename = splitPath[-1].split('.')
            PDBcode = splitFilename[-2].upper()
            fileName = PDBcode + '.' + splitFilename[-1]
            window_name = str(100*self.windowSize).replace('.', '_')
            pdt_name = str(self.PDT).replace('.', '_')
            PDBdirectory = 'Logfiles/%s_window_%s_pdt_%s/' % (
                PDBcode, window_name, pdt_name
                )
            file_name_start = '%s%s' % (PDBdirectory, PDBcode)
            newPathToInput = '%s%s' % (PDBdirectory, fileName)

            if os.path.isdir(PDBdirectory):
                print('Folder %s already exists locally at %s' % (
                    PDBdirectory, self.outputDir
                ))
                print('Do you want to overwrite the existing folder?\n'
                      '--USER INPUT-- type your choice and press RETURN\n'
                      'yes = overwrite this folder\n'
                      'no = do not overwrite this folder\n')
                owChoice = None
                while owChoice not in ['yes', 'no', 'y', 'n']:
                    if self.overwrite is True:
                        owChoice = 'yes'
                    elif self.overwrite is False and self.batchRun is True:
                        owChoice = 'no'
                    elif self.overwrite is False and self.batchRun is False:
                        owChoice = user_input(prompt).lower()

                    if owChoice == 'yes' or owChoice == 'y':
                        print('\nOverwriting existing folder')
                        shutil.rmtree(PDBdirectory)
                        copy_input(pathToInput, disk, newPathToInput,
                                   PDBdirectory)
                        break
                    elif owChoice == 'no' or owChoice == 'n':
                        print('\nKeeping original folder\n')
                        if self.batchRun is False:
                            sys.exit()
                        elif self.batchRun is True:
                            return
                        break
                    else:
                        print('Unrecognised input - please answer "yes" '
                              'or "no"')
            else:
                copy_input(pathToInput, disk, newPathToInput, PDBdirectory)

            # Checks that PDB / mmCIF file has been successfully copied to
            # the 'Logfiles' directory
            if not os.path.exists(newPathToInput):
                shutil.rmtree('%s' % PDBdirectory)
                print('\n\nERROR: Failed to copy %s to %s/%s/\n'
                      'Check that supplied file is not in use by another '
                      'program' % (pathToInput, self.outputDir, PDBdirectory))
                if self.batchRun is False:
                    sys.exit()
                elif self.batchRun is True:
                    return

            pathToInput = newPathToInput

        print('\nAll files generated by this program will be stored in:\n'
              '%s/%s/\n' % (self.outputDir, PDBdirectory))

        # Extracts ATOM/HETATM records from input PDB/mmCIF format file.
        exit = False
        if pathToInput[-4:] == '.cif':
            (atoms_list, disulfide_bonds, seqres, cryst1_line, exit
            ) = parse_mmcif_file(pathToInput)
        elif pathToInput[-4] == '.pdb':
            (atoms_list, disulfide_bonds, seqres, cryst1_line, exit
            ) = parse_pdb_file(pathToInput)

        if exit is True:
            shutil.rmtree('%s' % PDBdirectory)
            if self.batchRun is False:
                sys.exit()
            elif self.batchRun is True:
                return

        # Processes the ATOM/HETATM records to remove hydrogen atoms,
        # anisotropic B-factor records, and atoms with zero occupancy, as well
        # as retaining only the most probable alternate conformer for each atom.
        # These records are written to a PDB file that is fed into PDBCUR (a
        # program in the CCP4 softare suite); PDBCUR generates the unit cell
        # assembly from the coordinates of the processed asymmetric unit in the
        # input PDB file.
        print('\nProcessing input file to remove hydrogen atoms, anisotropic '
              '\nB factor records, and atoms with zero occupancy, as well as '
              '\nretaining only the most probable alternate conformations\n')
        exit = False
        pause = False
        exit, pause, clean_au_list, clean_au_file = clean_atom_rec(
            atoms_list, disulfide_bonds, seqres, cryst1_line, file_name_start
        )

        if exit is True:
            shutil.rmtree('%s' % PDBdirectory)
            if self.batchRun is False:
                sys.exit()
            elif self.batchRun is True:
                return

        if pause is True:
            print('\nTo enable damage detection, all macromolecule atoms '
                  'should be refined as full occupancy (i.e. 1),\n'
                  'except for alternate conformers (whose occupancy should sum '
                  'to 1).\nDisulfide bonds should be refined as single '
                  'conformers.\nOne or more atoms in the structure does not '
                  'meet these requirements (see ERRORs listed above).\n'
                  'Continue with RABDAM run?\n'
                  '--USER INPUT-- type your choice and press RETURN\n'
                  'yes = continue RABDAM run\n'
                  'no = terminate RABDAM run\n')
            owChoice = None
            while owChoice not in ['yes', 'no', 'y', 'n']:
                if self.batchRun is True:
                    owChoice = 'yes'  # Changed to yes for batch run on PDB only - change back to no before release!
                else:
                    owChoice = user_input(prompt).lower()

                if owChoice == 'yes' or owChoice == 'y':
                    print('Continuing RABDAM run\n')
                    break
                elif owChoice == 'no' or owChoice == 'n':
                    shutil.rmtree('%s' % PDBdirectory)
                    if self.batchRun is False:
                        sys.exit('\nExiting RABDAM\n')
                    elif self.batchRun is True:
                        return
                    break
                else:
                    print('Unrecognised input - please answer "yes" or "no"')

        # Deletes input .cif/.pdb file in Logfiles/ directory (i.e. the
        # downloaded/copied file, not the original file!) fed into the program,
        # unless createOrigpdb is set equal to True in the input file (default=False).
        if not self.createOrigpdb:
            os.remove(pathToInput)

        print('\nGenerating unit cell\n')
        PDBCURinputFile = '%sPDBCURinput.txt' % file_name_start
        PDBCURlog = '%sPDBCURlog.txt' % file_name_start
        genPDBCURinputs(PDBCURinputFile)
        unit_cell_pdb = '%s_unit_cell.pdb' % file_name_start
        runPDBCUR(clean_au_file, unit_cell_pdb, PDBCURinputFile, PDBCURlog)

        if not os.path.exists(unit_cell_pdb):
            print('\n\nERROR: Error in running PDBCUR, failed to generate unit '
                  'cell PDB file')
            shutil.rmtree('%s' % PDBdirectory)
            if self.batchRun is False:
                sys.exit()
            elif self.batchRun is True:
                return

        # The PDB file of the processed asymmetric unit is deleted unless
        # createAUpdb is set equal to True (default = False) in the input file.
        if self.createAUpdb is False:
            os.remove(clean_au_file)

        print('****************** End of Process PDB Section ******************\n'
              '****************************************************************\n')

        print('****************************************************************\n'
              '********************** Parsing PDB Section *********************\n')
        # Parses the newly generated unit cell PDB file into RABDAM, returning
        # a list of all atoms in the unit cell, plus the unit cell parameters.
        print('Reading in unit cell coordinates')
        ucAtomList = full_atom_list(unit_cell_pdb)

        # Halts program if no atoms selected for BDamage analysis
        if len(ucAtomList) < 1:
            print('\n\nERROR: No atoms selected for BDamage calculation')
            shutil.rmtree('%s' % PDBdirectory)
            if self.batchRun is False:
                sys.exit()
            elif self.batchRun is True:
                return

        # The unit cell PDB file is deleted unless createUCpdb is equal to True
        # in the input file (default = False).
        if self.createUCpdb is False:
            os.remove(unit_cell_pdb)

        print('\n****************** End of Parsing PDB Section ******************\n'
              '****************************************************************\n')

        print('****************************************************************\n'
              '****************** Translate Unit Cell Section *****************\n')
        # The unit cell parameters are converted into Cartesian vectors. These
        # vectors are then used to translate the unit cell -/+ 1 units in all
        # 3 (a, b and c) dimensions to generate a 3x3 parallelepiped. A PDB
        # file of this 3x3 parallelepiped is output if createAUCpdb is set
        # equal to True in the input file (default = False).

        transAtomList = np.empty([len(ucAtomList)*27, 3])
        transAtomIDList = ['']*(len(ucAtomList)*27)

        unit_cell_params = extract_unit_cell_params(cryst1_line)
        cartesianVectors = convertToCartesian(unit_cell_params)

        atom_count = 0
        for a in range(-1, 2):
            for b in range(-1, 2):
                for c in range(-1, 2):
                    transAtomList, transAtomIDList, atom_count = translateUnitCell(
                        ucAtomList, transAtomList, transAtomIDList,
                        cartesianVectors, a, b, c, atom_count,
                        self.createAUCpdb, self.createTApdb
                    )
        # Halts program if error in unit cell translation
        if atom_count != len(ucAtomList)*27:
            print('\n\nERROR: Failed to translate all unit cell atoms to '
                  'create 3x3 unit cell assembly')
            shutil.rmtree('%s' % PDBdirectory)
            if self.batchRun is False:
                sys.exit()
            elif self.batchRun is True:
                return

        # Creates PDB file of 3x3 unit cell assembly. WARNING: VERY slow and
        # RAM-consuming for large structures!
        if self.createAUCpdb is True:
            aucPDBfilepath = '%s_all_unit_cells.pdb' % file_name_start
            auc_pdb_atom_list = convert_array_to_atom_list(
                transAtomList, transAtomIDList, ucAtomList
            )
            makePDB('', auc_pdb_atom_list, '', seqres, aucPDBfilepath, 'Bfactor')

        print('\n************** End of Translate Unit Cell Section **************\n'
              '****************************************************************\n')

        print('****************************************************************\n'
              '********************* Trim Crystal Section *********************\n')
        # Atoms in the 3x3 parallelepiped which are a distance greater than the
        # packing density threshold from the boundaries of the processed
        # asymmetric unit are discarded. A PDB file of the trimmed
        # parallelepiped is created if createTApdb is set equal to True
        # (default = False) in the input file.
        bdamAtomList = b_damage_atom_list(
            clean_au_list, seqres, self.HETATM, self.protOrNA, self.addAtoms,
            self.removeAtoms
        )

        # Halts program if no atoms selected for BDamage analysis
        if len(bdamAtomList) < 1:
            print('\n\nERROR: No atoms selected for BDamage calculation')
            shutil.rmtree('%s' % PDBdirectory)
            if self.batchRun is False:
                sys.exit()
            elif self.batchRun is True:
                return

        auParams = getAUparams(bdamAtomList)
        print('\nObtained asymmetric unit parameters:')
        print('xMin = %8.3f' % auParams[0])
        print('xMax = %8.3f' % auParams[1])
        print('yMin = %8.3f' % auParams[2])
        print('yMax = %8.3f' % auParams[3])
        print('zMin = %8.3f' % auParams[4])
        print('zMax = %8.3f\n' % auParams[5])

        print('Removing atoms outside of packing density threshold')
        keepParams = convertParams(auParams, self.PDT)
        trimmedAtomList, trimmedAtomIDList = trimAtoms(
            transAtomList, keepParams, transAtomIDList, self.createAUCpdb,
            self.createTApdb, str(self.PDT)
        )

        # Creates PDB file of trimmed 3x3 unit cell assembly. WARNING: VERY
        # slow and RAM-consuming for large structures!
        if self.createTApdb is True:
            taPDBfilepath = '%s_trimmed_atoms.pdb' % file_name_start
            ta_pdb_atom_list = convert_array_to_atom_list(
                trimmedAtomList, trimmedAtomIDList, ucAtomList
            )
            makePDB('', ta_pdb_atom_list, '', seqres, taPDBfilepath, 'Bfactor')

        print('\n****************** End of Trim Crystal Section *****************\n'
              '****************************************************************\n')

        print('****************************************************************\n'
              '*************** Calculate Packing Density Section **************\n')
        # Calculates the packing density (atomic contact number) of every atom
        # in the asymmetric unit.

        print('Calculating packing density values\n')
        au_atom_xyz = get_xyz_from_objects(bdamAtomList)
        packing_density_array = calc_packing_density(
            au_atom_xyz, trimmedAtomList, self.PDT
        )
        write_pckg_dens_to_atoms(bdamAtomList, packing_density_array)

        print('*********** End of Calculate Packing Density Section ***********\n'
              '****************************************************************\n')

        print('****************************************************************\n'
              '****************** Calculate BDamage Section *******************\n')
        # Atoms in the asymmetric unit are ordered via their packing density
        # values; the BDamage value of each atom is then calculated as the
        # ratio of its B-factor as compared to the average of the B-factor
        # values of similarly (identified via sliding window) packed atoms.

        print('Calculating BDamage values\n')
        window = int(round((len(bdamAtomList)*self.windowSize), 0))
        if window % 2 == 0:
            window = window + 1  # Window size must be an odd number.
        if window < 11:
            window = 11  # Minimum window size is 11.
        print('Size of sliding window --> %s atoms\n' % window)
        calcBDam(bdamAtomList, window)

        print('****************************************************************\n'
              '******************* Writing DataFrame Section ******************\n')
        # Writes properties of atoms to be considered for BDamage analysis to a
        # DataFrame. The DataFrame, plus additional variables and lists
        # required for subsequent analysis, are then pickled - this allows
        # multiple analysis runs to be performed from the same DataFrame,
        # thereby reducing their calculation time.

        print('Writing BDamage data to DataFrame\n')
        df = writeDataFrame(bdamAtomList)
        print('Saving DataFrame\n')
        storage = '%s/DataFrame' % PDBdirectory
        os.mkdir(storage)
        storage_file = '%s/%s' % (storage, PDBcode)
        df.to_pickle(storage_file + '_dataframe.pkl')
        with open(storage_file + '_variables.pkl', 'wb') as f:
            pickle.dump((file_name_start, seqres, bdamAtomList, window), f)

        print('****************************************************************\n'
              '*************** End Of Writing DataFrame Section ***************\n')

        # Changes directory back to the 'RABDAM' directory (that in which the
        # rabdam.py script is saved).
        os.chdir('%s' % cwd)

    def rabdam_analysis(self, run, output_options):
        """
        Uses values in DataFrame returned from calling the 'rabdam_dataframe'
        function to write output analysis files.
        """

        prompt = '> '
        import os
        import sys
        import pickle
        import pandas as pd

        if sys.version_info[0] < 3:
            user_input = raw_input
        else:
            user_input = input

        if __name__ == 'Subroutines.CalculateBDamage':
            from Subroutines.output import generate_output_files
            from Subroutines.makeDataFrame import makePDB
        else:
            from rabdam.Subroutines.output import generate_output_files
            from rabdam.Subroutines.makeDataFrame import makePDB

        if run == 'rabdam_analysis':
            print('************************ RABDAM ANALYSIS ***********************\n')

        # Changes directory to the specified location for the output 'Logfiles'
        # directory. The default location is the current working directory
        # (i.e. that in which the rabdam.py script is saved).
        cwd = os.getcwd()
        os.chdir('%s' % self.outputDir)

        print('\n****************************************************************\n'
              '***************** Processing DataFrame Section *****************\n')
        # Checks that pkl files output by 'rabdam_dataframe' function exist.
        # Then checks if output directory specified in input file (default =
        # current working directory) already contains any analysis output files
        # - if it does then user input is requested ('Do you want to overwrite
        # the existing analysis files?'):
        # yes = all old analysis files are removed and replaced by new analysis
        #       files
        # no = old analysis files are retained, exit program
        # Note that currently there is no option to replace only a subset of
        # the output analysis files.

        pathToInput = self.pathToInput.replace('\\', '/')
        splitPath = pathToInput.split('/')
        PDBcode = splitPath[len(splitPath)-1]
        for file_type in ['.pdb', '.cif']:
            PDBcode = PDBcode.replace('%s' % file_type, '')
        PDBcode = PDBcode.upper()
        window_name = str(100*self.windowSize).replace('.', '_')
        pdt_name = str(self.PDT).replace('.', '_')
        PDBdirectory = 'Logfiles/%s_window_%s_pdt_%s' % (PDBcode, window_name,
                                                         pdt_name)
        PDB_analysis_file = '%s/%s' % (PDBdirectory, PDBcode)
        storage_directory = '%s/DataFrame' % PDBdirectory
        storage_file = '%s/%s' % (storage_directory, PDBcode)

        if not os.path.isdir(storage_directory):
            print('\n\nERROR: Folder %s not found\n' % storage_directory)
            if self.batchRun is False:
                sys.exit('Exiting RABDAM\n')
            elif self.batchRun is True:
                return

        potential_analysis_files = [
            '_BDamage.csv', '_BDamage.pdb', '_BDamage.cif', '_BDamage.svg',
            '_Bnet_protein.svg', '_Bnet_protein.pkl', '_Bnet_NA.svg',
            '_Bnet_NA.pkl', '_BDamage.html'
        ]
        actual_analysis_files = []
        for file_name in potential_analysis_files:
            if os.path.isfile(PDB_analysis_file + file_name):
                actual_analysis_files.append(PDB_analysis_file + file_name)
        if len(actual_analysis_files) > 0:
            print('There are one or more RABDAM analysis files already present\n'
                  'in folder %s' % PDBdirectory)
            print('Do you want to overwrite the existing analysis files?\n'
                  '--USER INPUT-- type your choice and press RETURN\n'
                  'yes = overwrite ALL analysis files\n'
                  'no = do not overwrite analysis files\n')
            owChoice = None
            while owChoice not in ['yes', 'no', 'y', 'n']:
                if self.overwrite is True:
                    owChoice = 'yes'
                elif self.overwrite is False and self.batchRun is True:
                    owChoice = 'no'
                elif self.overwrite is False and self.batchRun is False:
                    owChoice = user_input(prompt).lower()

                if owChoice == 'yes' or owChoice == 'y':
                    print('\nOverwriting existing analysis files\n')
                    for file_name in actual_analysis_files:
                        os.remove(file_name)
                    break
                elif owChoice == 'no' or owChoice == 'n':
                    print('Keeping original analysis files\nExiting RABDAM\n')
                    if self.batchRun is False:
                        sys.exit()
                    elif self.batchRun is True:
                        return
                    break
                else:
                    print('Unrecognised input - please answer "yes" or "no"')

        # Pkl files unpickled
        df = pd.read_pickle(storage_file + '_dataframe.pkl')
        print('Unpickling DataFrame and variables\n')
        with open(storage_file + '_variables.pkl', 'rb') as f:
            (file_name_start, seqres, bdamAtomList, window) = pickle.load(f)

        print('************** End Of Processing DataFrame Section *************\n'
              '****************************************************************\n')

        print('****************************************************************\n'
              '***************** Writing Output Files Section *****************\n')
        # Uses the values stored in the BDamage DataFrame to:
        # - write DataFrame values to a csv file (providing users with a copy
        #   of the raw data that they can manipulate for further analysis as
        #   they wish)
        # - write a PDB file in which the B-factors are replaced by BDamage
        #   values (allowing users to e.g. colour the structure by BDamage
        #   when viewed with molecular graphics software)
        # - write a cif file in which a column listing the calculated BDamage
        #   values has been appended to the ATOM (/HETATM) records
        # - plot a kernel density estimate of the complete (i.e. all analysed
        #   atoms) BDamage distribution
        # - plot a kernel density estimate of the BDamage values of the
        #   carboxyl group oxygens of Asp and Glu residues, this plot is then
        #   used to calculate the value of the Bnet summary metric

        output = generate_output_files(
            out_file_start=file_name_start, pdb_code=PDBcode, df=df
        )

        if 'csv' in output_options:
            print('\nWriting csv file')
            output.make_csv(window)

        if 'pdb' in output_options:
            print('\nWriting PDB file with BDamage values replacing Bfactors')
            pdb_file_name = file_name_start + '_BDamage.pdb'
            makePDB('', bdamAtomList, '', seqres, pdb_file_name, 'BDamage')

        if 'cif' in output_options:
            print('\nWriting cif file with BDamage column')
            output.write_output_cif(bdamAtomList)

        if 'kde' in output_options or 'summary' in output_options:
            print('\nPlotting kernel density estimate')
            output.make_histogram(self.highlightAtoms)

        if 'bnet' in output_options or 'summary' in output_options:
            print('\nCalculating Bnet')
            output.calculate_Bnet(window_name, pdt_name, window)

        if 'summary' in output_options:
            print('\nWriting summary html file\n')
            output.write_html_summary(output_options, self.highlightAtoms)

        print('************** End of Writing Output Files Section *************\n'
              '****************************************************************\n')

        # Changes directory back to the 'RABDAM' directory (i.e. that in which
        # the rabdam.py script is saved).
        os.chdir('%s' % cwd)
