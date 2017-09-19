
# RABDAM
# Copyright (C) 2017 Garman Group, University of Oxford

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


class rabdam():
    def __init__(self, pathToPDB, outputDir, PDT, windowSize, protOrNA, HETATM,
                 addAtoms, removeAtoms, highlightAtoms, createOrigpdb,
                 createAUpdb, createUCpdb, createAUCpdb,
                 createTApdb):
        self.pathToPDB = pathToPDB
        self.outputDir = outputDir
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
        # Calculates BDamage for selected atoms within input PDB file and
        # writes output to DataFrame.

        prompt = '> '
        import sys
        import os
        import shutil
        import copy
        duplicate = copy.copy
        import pickle

        from PDBCUR import clean_pdb_file, genPDBCURinputs, runPDBCUR
        from parsePDB import (full_atom_list, b_damage_atom_list, downloadPDB,
                              copyPDB)
        from translateUnitCell import convertToCartesian, translateUnitCell
        from trimUnitCellAssembly import getAUparams, convertParams, trimAtoms
        from makeDataFrame import makePDB, writeDataFrame
        from Bdamage import (calcBdam, get_xyz_from_objects,
                             calc_packing_density, write_pckg_dens_to_atoms)

        if run == 'rabdam':
            print '**************************** RABDAM ****************************\n'
            print('\nPlease cite: M. Gerstel, C. M. Deane and E.F. Garman. (2015).\n'
                  'J. Synchrotron Radiation, 22, 201-212.\n'
                  'http://dx.doi.org/doi:10.1107/S1600577515002131\n')
        elif run == 'rabdam_dataframe':
            print '*********************** RABDAM DATAFRAME ***********************\n'
            print('\nPlease cite: M. Gerstel, C. M. Deane and E.F. Garman. (2015).\n'
                  'J. Synchrotron Radiation, 22, 201-212.\n'
                  'http://dx.doi.org/doi:10.1107/S1600577515002131\n')

        print('\n****************************************************************\n'
              '***************** Program to calculate BDamage *****************\n'
              '****************************************************************\n')

        print('****************************************************************\n'
              '************************* Input Section ************************\n')
        # Prints the values of the program options to be used in the current
        # RABDAM run
        print 'Calculating BDamage for %s' % self.pathToPDB
        print 'Writing output files to %s' % self.outputDir
        if self.PDT == 7:
            print 'Using default packing density threshold of 7 Angstroms'
        else:
            print 'Packing density threshold defined by user as %s Angstroms' % self.PDT
        if self.windowSize == 0.02:
            print 'Using default window size of 2%'
        else:
            print 'Window size defined by user as %s%%' % (self.windowSize*100)
        if self.HETATM is True:
            print 'Keeping HETATM'
        elif self.HETATM is False:
            print 'Removing HETATM'
        if self.protOrNA == 'PROTEIN':
            print 'Retaining protein atoms, discarding nucleic acid atoms'
        elif self.protOrNA in ['NUCLEICACID', 'NA']:
            print 'Retaining nucleic acid atoms, discarding protein atoms'
        if len(self.addAtoms) == 0:
            print 'No atoms to be added'
        else:
            add_atoms_string = ''
            for value in self.addAtoms:
                value = value + ', '
                add_atoms_string = add_atoms_string + value
            add_atoms_string = add_atoms_string.strip(', ')
            print 'Atoms to be added: %s' % add_atoms_string
        if len(self.removeAtoms) == 0:
            print 'No atoms to be removed'
        else:
            remove_atoms_string = ''
            for value in self.removeAtoms:
                value = value + ', '
                remove_atoms_string = remove_atoms_string + value
            remove_atoms_string = remove_atoms_string.strip(', ')
            print 'Atoms to be removed: %s' % remove_atoms_string

        print('\n********************* End of Input Section *********************\n'
              '****************************************************************\n')

        # Changes directory to the specified location for the output 'Logfiles'
        # directory. The default location is the current working directory
        # (i.e. that in which the rabdam.py script is saved).
        cwd = os.getcwd()
        os.chdir('%s' % self.outputDir)

        print('****************************************************************\n'
              '********************** Process PDB Section *********************\n')
        # Creates a new directory named after the input PDB file (after
        # checking that this directory does not already exist) in the
        # 'Logfiles' directory. Then saves a copy of the input PDB file to the
        # new directory.

        # If 4 digit PDB accession code has been supplied:
        if len(self.pathToPDB) == 4:
            print 'PDB code supplied\n'
            PDBcode = self.pathToPDB.upper()
            window_name = 100*self.windowSize
            window_name = str(window_name).replace('.', '_')
            pdt_name = str(self.PDT).replace('.', '_')
            PDBdirectory = 'Logfiles/%s_window_%s_pdt_%s/' % (PDBcode, window_name, pdt_name)
            pdb_file_path = '%s%s' % (PDBdirectory, PDBcode)
            pathToPDB = '%s%s.pdb' % (PDBdirectory, PDBcode)

            # If directory with same name as PDBdirectory already exists in
            # 'Logfiles' directory, user input is requested ('Do you want to
            # overwrite the existing folder?'):
            # yes = old PDBdirectory is deleted, new PDBdirectory is created
            #       and copy of the PDB file is downloaded from the RSCB PDB
            #       website and saved to the new directory
            # no = old PDBdirectory is retained, exit program
            # If it doesn't already exist, new PDBdirectory is created and
            # copy of the PDB file is downloaded from the RSCB PDB website and
            # saved to the new directory.
            if os.path.isdir(PDBdirectory):
                print 'Folder %s already exists locally at %s' % (PDBcode,
                                                                  PDBdirectory)
                print('Do you want to overwrite the existing folder?\n'
                      '--USER INPUT-- type your choice and press RETURN\n'
                      'yes = overwrite this folder\n'
                      'no = do not overwrite this folder\n')
                owChoice = None
                while owChoice not in ['yes', 'no', 'y', 'n']:
                    owChoice = raw_input(prompt).lower()
                    if owChoice == 'yes' or owChoice == 'y':
                        print '\nOverwriting existing folder'
                        shutil.rmtree(PDBdirectory)
                        downloadPDB(PDBcode, PDBdirectory, pathToPDB)
                        break
                    elif owChoice == 'no' or owChoice == 'n':
                        print('\nKeeping original folder\n'
                              'Exiting RABDAM')
                        sys.exit()
                        break
                    else:
                        print 'Unrecognised input - please answer "yes" or "no"'
            else:
                downloadPDB(PDBcode, PDBdirectory, pathToPDB)

            # Checks that PDB file has been successfully downloaded and saved to
            # the 'Logfiles' directory
            if not os.path.exists(pathToPDB):
                sys.exit('ERROR: Failed to download and save PDB file - cause unknown')

        # If filepath to PDB has been supplied:
        else:
            # Changes directory to allow input PDB file (a '.pdb' or '.txt'
            # file) to be read from any provided file path. If directory with
            # same name as PDBdirectory already exists in 'Logfiles' directory,
            # user input is requested ('Do you want to overwrite the existing
            # folder?'):
            # yes = old PDBdirectory is deleted, new PDBdirectory is created
            #       and copy of input PDB file is saved to the new directory
            # no = old PDBdirectory is retained, exit program
            # If it doesn't already exist, new PDBdirectory is created and copy
            # of input PDB file is saved to the new directory.
            owd = os.getcwd()
            pathToPDB = self.pathToPDB.replace('\\', '/')
            splitPath = pathToPDB.split('/')
            disk = '%s/' % splitPath[0]
            os.chdir('/')
            os.chdir(disk)
            if not os.path.exists(pathToPDB):
                sys.exit('ERROR: Supplied filepath not recognised')
            os.chdir(owd)

            if pathToPDB[-4:] not in ['.pdb', '.txt']:
                sys.exit('ERROR: Supplied filepath to PDB is not a .pdb or'
                         '.txt file')
            else:
                print 'Filepath to .pdb or .txt file supplied\n'
                fileName = splitPath[len(splitPath)-1]
                splitFilename = fileName.split('.')
                PDBcode = splitFilename[len(splitFilename)-2].upper()
                fileName = PDBcode + '.' + splitFilename[len(splitFilename)-1]
                window_name = 100*self.windowSize
                window_name = str(window_name).replace('.', '_')
                pdt_name = str(self.PDT).replace('.', '_')
                PDBdirectory = 'Logfiles/%s_window_%s_pdt_%s/' % (PDBcode, window_name, pdt_name)
                pdb_file_path = '%s%s' % (PDBdirectory, PDBcode)
                newPathToPDB = '%s%s' % (PDBdirectory, fileName)

                if os.path.isdir(PDBdirectory):
                    print 'Folder %s already exists locally at %s' % (
                        PDBcode, PDBdirectory)
                    print('Do you want to overwrite the existing folder?\n'
                          '--USER INPUT-- type your choice and press RETURN\n'
                          'yes = overwrite this folder\n'
                          'no = do not overwrite this folder\n')
                    owChoice = None
                    while owChoice not in ['yes', 'no', 'y', 'n']:
                        owChoice = raw_input(prompt).lower()
                        if owChoice == 'yes' or owChoice == 'y':
                            print '\nOverwriting existing folder'
                            shutil.rmtree(PDBdirectory)
                            copyPDB(pathToPDB, disk, newPathToPDB, PDBdirectory)
                            break
                        elif owChoice == 'no' or owChoice == 'n':
                            print('\nKeeping original folder\n'
                                  'Exiting RABDAM')
                            sys.exit()
                            break
                        else:
                            print 'Unrecognised input - please answer "yes" or "no"'
                else:
                    copyPDB(pathToPDB, disk, newPathToPDB, PDBdirectory)

                # Checks that PDB file has been successfully copied to the
                # 'Logfiles' directory
                if not os.path.exists(newPathToPDB):
                    sys.exit('ERROR: Failed to copy PDB file to Logfiles '
                             'directory.\nCheck that supplied PDB file is not '
                             'in use by another program')

                pathToPDB = newPathToPDB

        print('All files generated by this program will be stored in:\n'
              '    %s\n' % PDBdirectory)

        # Processes the input PDB file to remove hydrogen atoms, anisotropic
        # B factor records, and atoms with zero occupancy, as well as
        # retaining only the most probable alternate conformations. The unit
        # cell assembly is then generated from the coordinates of the processed
        # PDB file and its associated symmetry operations. Note that unit cell
        # generation is currently performed by the PDBCUR program from the CCP4
        # software suite.
        print ('\nProcessing PDB file to remove hydrogen atoms, anisotropic '
               '\nB factor records, and atoms with zero occupancy, as well as '
               '\nretaining only the most probable alternate conformations')

        (clean_au_file, clean_au_list, header_lines, footer_lines,
         unit_cell_params) = clean_pdb_file(pathToPDB, pdb_file_path)

        # Deletes input PDB file fed into the program unless createOrigpdb
        # is set equal to True in the input file (default=False).
        if not self.createOrigpdb:
            os.remove(pathToPDB)

        print '\nGenerating unit cell\n'
        PDBCURinputFile = '%sPDBCURinput.txt' % pdb_file_path
        PDBCURlog = '%sPDBCURlog.txt' % pdb_file_path
        genPDBCURinputs(PDBCURinputFile)
        unit_cell_pdb = '%s_unit_cell.pdb' % pdb_file_path
        runPDBCUR(clean_au_file, unit_cell_pdb, PDBCURinputFile, PDBCURlog)

        if not os.path.exists(unit_cell_pdb):
            sys.exit('ERROR: Error in running PDBCUR, failed to generate Unit'
                     'Cell PDB file')

        print('****************** End of Process PDB Section ******************\n'
              '****************************************************************\n')

        print('****************************************************************\n'
              '********************** Parsing PDB Section *********************\n')
        # Parses the newly generated unit cell PDB file into RABDAM, returning
        # a list of all atoms in the unit cell, plus the unit cell parameters.
        # The unit cell PDB file is then deleted unless createUCpdb is set
        # equal to True in the input file (default = False).
        print 'Reading in unit cell coordinates'
        ucAtomList = full_atom_list(unit_cell_pdb)

        if self.createUCpdb is False:
            os.remove(unit_cell_pdb)

        print('\n****************** End of Parsing PDB Section ******************\n'
              '****************************************************************\n')

        print('****************************************************************\n'
              '****************** Translate Unit Cell Section *****************\n')
        # The unit cell parameters are converted into Cartesian vectors. These
        # vectors are then used to translate the unit cell -/+ 1 units in all
        # (3 - a, b and c) dimensions to generate a 3x3 parallelepiped. A PDB
        # file of this 3x3 parallelepiped is output if createAUCpdb is set
        # equal to True in the input file (default = False).

        transAtomList = duplicate(ucAtomList)
        cartesianVectors = convertToCartesian(unit_cell_params)
        for a in xrange(-1, 2):
            for b in xrange(-1, 2):
                for c in xrange(-1, 2):
                    if a == 0 and b == 0 and c == 0:
                        pass  # No identity translation
                    else:
                        transAtomList = translateUnitCell(ucAtomList,
                                                          transAtomList,
                                                          cartesianVectors,
                                                          a, b, c)

        if self.createAUCpdb is True:
            aucPDBfilepath = '%s_all_unit_cells.pdb' % pdb_file_path
            makePDB(header_lines, transAtomList, footer_lines, aucPDBfilepath,
                    'Bfactor')

        print('\n************** End of Translate Unit Cell Section **************\n'
              '****************************************************************\n')

        print('****************************************************************\n'
              '********************* Trim Crystal Section *********************\n')
        # Atoms in the 3x3 parallelepiped which are a distance greater than the
        # packing density threshold from the boundaries of the processed
        # asymmetric unit are discarded. A PDB file of the trimmed
        # parallelepiped is created if createTApdb is set equal to True
        # (default = False) in the input file. The PDB file of the processed
        # asymmetric unit is deleted unless createAUpdb is set equal to True
        # (default = False) in the input file.
        bdamAtomList = b_damage_atom_list(clean_au_list, self.HETATM,
                                          self.protOrNA, self.addAtoms,
                                          self.removeAtoms)

        if self.createAUpdb is False:
            os.remove(clean_au_file)

        auParams = getAUparams(bdamAtomList)
        print '\nObtained asymmetric unit parameters:'
        print 'xMin = %8.3f' % auParams[0]
        print 'xMax = %8.3f' % auParams[1]
        print 'yMin = %8.3f' % auParams[2]
        print 'yMax = %8.3f' % auParams[3]
        print 'zMin = %8.3f' % auParams[4]
        print 'zMax = %8.3f\n' % auParams[5]

        print 'Removing atoms outside of packing density threshold'
        keepParams = convertParams(auParams, self.PDT)
        trimmedAtomList = trimAtoms(transAtomList, keepParams)

        if self.createTApdb is True:
            taPDBfilepath = '%s_trimmed_atoms.pdb' % pdb_file_path
            makePDB(header_lines, trimmedAtomList, footer_lines, taPDBfilepath,
                    'Bfactor')

        print('\n****************** End of Trim Crystal Section *****************\n'
              '****************************************************************\n')

        print('****************************************************************\n'
              '*************** Calculate Packing Density Section **************\n')
        # Calculates the packing density (atomic contact number) of every atom
        # in the asymmetric unit.

        print 'Calculating packing density values\n'
        au_atom_xyz, trim_atom_xyz = get_xyz_from_objects(bdamAtomList,
                                                          trimmedAtomList)
        packing_density_array = calc_packing_density(au_atom_xyz,
                                                     trim_atom_xyz, self.PDT)
        write_pckg_dens_to_atoms(bdamAtomList, packing_density_array)

        print('*********** End of Calculate Packing Density Section ***********\n'
              '****************************************************************\n')

        print('****************************************************************\n'
              '****************** Calculate BDamage Section *******************\n')
        # Atoms in the asymmetric unit are ordered via their packing density
        # values; the BDamage value of each atom is then calculated as the
        # ratio of its Bfactor as compared to the average of the Bfactor values
        # of similarly (identified via sliding window) packed atoms.

        print 'Calculating BDamage values\n'
        window = int(round((len(bdamAtomList)*self.windowSize), 0))
        if window % 2 == 0:
            window = window + 1  # Window size must be an odd number.
        if window < 3:
            window == 3  # Minimum window size is 3.
        print 'Size of sliding window --> %s atoms\n' % window
        calcBdam(bdamAtomList, window)

        print('****************************************************************\n'
              '******************* Writing DataFrame Section ******************\n')
        # Writes properties of atoms to be considered for BDamage analysis to a
        # DataFrame. The DataFrame, plus additional variables and lists
        # required for subsequent analysis, are then pickled - this allows
        # multiple analysis runs to be performed from the same DataFrame,
        # thereby reducing their calculation time.

        print 'Writing BDamage data to DataFrame\n'
        df = writeDataFrame(bdamAtomList)

        print 'Saving DataFrame\n'
        storage = '%s/DataFrame' % PDBdirectory
        os.makedirs(storage)
        storage_file = '%s/%s' % (storage, PDBcode)
        df.to_pickle(storage_file + '_dataframe.pkl')
        with open(storage_file + '_variables.pkl', 'wb') as f:
            pickle.dump((pdb_file_path, PDBcode, bdamAtomList, header_lines,
                         footer_lines, window), f)

        print('****************************************************************\n'
              '*************** End Of Writing DataFrame Section ***************\n')

        # Changes directory back to the 'RABDAM' directory (that in which the
        # rabdam.py script is saved).
        os.chdir('%s' % cwd)

    def rabdam_analysis(self, run, output_options, count):
        # Uses values in DataFrame returned from calling the 'rabdam_dataframe'
        # function to write output analysis files.

        import os
        import sys
        prompt = '> '
        import pickle
        import pandas as pd
        from output import generate_output_files
        from makeDataFrame import makePDB

        if run == 'rabdam_analysis':
            print '************************ RABDAM ANALYSIS ***********************\n'
            print('\nPlease cite: M. Gerstel, C. M. Deane and E.F. Garman. (2015).\n'
                  'J. Synchrotron Radiation, 22, 201-212.\n'
                  'http://dx.doi.org/doi:10.1107/S1600577515002131\n')

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

        pathToPDB = self.pathToPDB.replace('\\', '/')
        splitPath = pathToPDB.split('/')
        pathToPDB = splitPath[len(splitPath)-1]
        PDBcode = pathToPDB.replace('.pdb', '')
        PDBcode = PDBcode.replace('.txt', '')
        PDBcode = PDBcode.upper()
        window_name = 100*self.windowSize
        window_name = str(window_name).replace('.', '_')
        pdt_name = str(self.PDT).replace('.', '_')
        PDBdirectory = 'Logfiles/%s_window_%s_pdt_%s' % (PDBcode, window_name,
                                                         pdt_name)
        PDB_analysis_file = '%s/%s' % (PDBdirectory, PDBcode)
        storage_directory = '%s/DataFrame' % PDBdirectory
        storage_file = '%s/%s' % (storage_directory, PDBcode)

        if not os.path.isdir(storage_directory):
            print 'Folder %s not found' % (storage_directory)
            print 'Exiting RABDAM analysis'
            sys.exit()

        potential_analysis_files = ['_BDamage.csv', '_BDamage.html',
                                    '_BDamage.pdb', '_BDamage.svg',
                                    '_Bnet_Protein.svg', '_Bnet_NA.svg']
        actual_analysis_files = []
        for name in potential_analysis_files:
            if os.path.isfile(PDB_analysis_file + name):
                actual_analysis_files.append(PDB_analysis_file + name)
        if len(actual_analysis_files) > 0:
            print('There are one or more RABDAM analysis files already present\n'
                  'in folder %s' % PDBdirectory)
            print('Do you want to overwrite the existing analysis files?\n'
                  '--USER INPUT-- type your choice and press RETURN\n'
                  'yes = overwrite ALL analysis files\n'
                  'no = do not overwrite analysis files\n')
            owChoice = None
            while owChoice not in ['yes', 'no', 'y', 'n']:
                owChoice = raw_input(prompt).lower()
                if owChoice == 'y' or owChoice == 'y':
                    print '\nOverwriting existing analysis files\n'
                    for name in actual_analysis_files:
                        os.remove(name)
                    break
                elif owChoice == 'no' or owChoice == 'n':
                    print('Keeping original analysis files\n'
                          'Exiting RABDAM')
                    sys.exit()
                    break
                else:
                    print 'Unrecognised input - please answer "yes" or "no"'

        # Pkl files unpickled
        print 'Unpickling DataFrame and variables\n'
        with open(storage_file + '_variables.pkl', 'rb') as f:
            (pdb_file_path, PDBcode, bdamAtomList, header_lines, footer_lines,
             window) = pickle.load(f)
        df = pd.read_pickle(storage_file + '_dataframe.pkl')

        print('************** End Of Processing DataFrame Section *************\n'
              '****************************************************************\n')

        print('****************************************************************\n'
              '***************** Writing Output Files Section *****************\n')
        # Uses the values stored in the BDamage DataFrame to:
        # - write DataFrame values to csv file (providing users with a copy of
        #   the raw data that they can manipulate for further analysis as they
        #   wish)
        # - write a PDB file in which the Bfactors are replaced by BDamage
        #   values (allowing users to e.g. colour the structure by BDamage
        #   when viewed with molecular graphics software)
        # - plot a kernel density estimate of the complete (i.e. all analysed
        #   atoms) BDamage distribution
        # - plot a kernel density estimate of the BDamage values of the
        #   carboxyl group oxygens of Asp and Glu residues, this plot is then
        #   used to calculate the value of the Bnet summary metric

        output = generate_output_files(pdb_file_path=pdb_file_path, df=df)

        if 'csv' in output_options:
            print 'Writing csv file\n'
            output.make_csv(bdamAtomList, window)

        if 'pdb' in output_options:
            print 'Writing PDB file with BDamage values replacing Bfactors'
            pdb_file_name = pdb_file_path + '_BDamage.pdb'
            makePDB(header_lines, bdamAtomList, footer_lines, pdb_file_name,
                    'BDamage')

        if 'kde' in output_options:
            print '\nPlotting kernel density estimate\n'
            output.make_histogram(self.highlightAtoms)

        if 'bnet' in output_options:
            print 'Calculating Bnet\n'
            output.calculate_Bnet(window_name, pdt_name, count, window)

        if 'summary' in output_options and 'kde' in output_options and 'bnet' in output_options:
            print 'Writing summary html file\n'
            output.write_html_summary(cwd, self.highlightAtoms)

        print('************** End of Writing Output Files Section *************\n'
              '****************************************************************\n')

        # Changes directory back to the 'RABDAM' directory (i.e. that in which
        # the rabdam.py script is saved).
        os.chdir('%s' % cwd)
