
class rabdam():
    def __init__(self, pathToPDB, outputDir, PDT, windowSize, protOrNA, HETATM,
                 addAtoms, removeAtoms, threshold, highlightAtoms, createAUpdb,
                 createUCpdb, createAUCpdb, createTApdb):
        self.pathToPDB = pathToPDB
        self.outputDir = outputDir
        self.PDT = PDT
        self.windowSize = windowSize
        self.protOrNA = protOrNA
        self.HETATM = HETATM
        self.addAtoms = addAtoms
        self.removeAtoms = removeAtoms
        self.threshold = threshold
        self.highlightAtoms = highlightAtoms
        self.createAUpdb = createAUpdb
        self.createUCpdb = createUCpdb
        self.createAUCpdb = createAUCpdb
        self.createTApdb = createTApdb

    def rabdam_dataframe(self, run):
        # Calculates B_Damage values from input PDB file and writes
        # output to DataFrame.

        prompt = '> '
        import sys
        import os
        import shutil
        import copy
        duplicate = copy.copy
        import pickle

        from PDBCUR import clean_pdb_file, genPDBCURinputs, runPDBCUR
        from parsePDB import (full_atom_list, b_damage_atom_list, downloadPDB,
                              copyPDB, getAUparams, trimAtoms)
        from translateUnitCell import (convertToCartesian, translateUnitCell)
        from makePDB import makePDB, writeBdam
        from atomCheck import convertParams
        from Bdamage import (calcBdam, get_xyz_from_objects,
                             calc_packing_density, write_pckg_dens_to_atoms)

        if run == 'rabdam':
            print '**************************** RABDAM ****************************\n'
            print('\nPlease cite: M. Gerstel, C. M. Deane and E.F. Garman. (2015).\n'
                  'J. Synchrotron Radiation, 22, 201-212\n'
                  'http://dx.doi.org/doi:10.1107/S1600577515002131\n')
        elif run == 'rabdam_dataframe':
            print '*********************** RABDAM DATAFRAME ***********************\n'
            print('\nPlease cite: M. Gerstel, C. M. Deane and E.F. Garman. (2015).\n'
                  'J. Synchrotron Radiation, 22, 201-212\n'
                  'http://dx.doi.org/doi:10.1107/S1600577515002131\n')

        print('\n****************************************************************\n'
              '**************** Program to calculate B_Damage *****************\n'
              '****************************************************************\n')

        print('****************************************************************\n'
              '************************* Input Section ************************\n')
        # Prints the input values read into the program from INPUT.txt.

        print 'Calculating B_Damage for %s' % self.pathToPDB
        print 'Writing output files to %s' % self.outputDir
        if self.PDT == 14:
            print 'Using default packing density threshold of 14 Angstroms'
        else:
            print 'Packing density threshold defined by user as %s Angstroms' % self.PDT
        if self.windowSize == 0.02:
            print 'Using default window size of 0.02'
        else:
            print 'Window size defined by user as %s' % self.windowSize
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
            for value in self.addAtoms:
                print 'Atoms to be added: %s' % value
        if len(self.removeAtoms) == 0:
            print 'No atoms to be removed'
        else:
            for value in self.removeAtoms:
                print 'Atoms to be removed: %s' % value

        print('\n********************* End of Input Section *********************\n'
              '****************************************************************\n')

        cwd = os.getcwd()
        os.chdir('%s' % self.outputDir)

        print('****************************************************************\n'
              '********************** Process PDB Section *********************\n')
        # Creates a new directory named after the input PDB file (after checking
        # that this directory does not already exist) in the Logfiles directory.
        # Then saves a copy of the input PDB file to the new directory.

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

            # If PDB directory of this name already exists in Logfiles directory,
            # user input requested ('Do you want to overwrite the existing file?'):
            # yes = old PDB directory is deleted, new PDB directory is created and
            #       copy of the PDB file is downloaded from the RSCB PDB website
            #       and saved to the new directory
            # no = old PDB directory is retained, exit program
            # If it doesn't already exist, new PDB directory is created and copy of
            # the PDB file is downloaded from the RSCB PDB website and saved to
            # the new directory.
            if os.path.isdir(PDBdirectory):
                print 'Folder %s already exists locally at %s' % (PDBcode,
                                                                  PDBdirectory)
                print('Do you want to overwrite the existing file?\n'
                      '--USER INPUT-- type your choice and press RETURN\n'
                      'yes = overwrite this folder\n'
                      'no = do not overwrite this folder\n')
                owChoice = None
                while owChoice not in ['YES', 'NO', 'Y', 'N']:
                    owChoice = raw_input(prompt).upper()
                    if owChoice == 'YES' or owChoice == 'Y':
                        print '\nOverwriting existing folder'
                        shutil.rmtree(str(PDBdirectory))
                        downloadPDB(PDBcode, PDBdirectory, pathToPDB)
                        break
                    elif owChoice == 'NO' or owChoice == 'N':
                        print('\nKeeping original folder\n'
                              'Exiting RABDAM')
                        sys.exit()
                        break
                    else:
                        print 'Unrecognised input - please answer "yes" or "no"'
            else:
                downloadPDB(PDBcode, PDBdirectory, pathToPDB)

            # Checks that PDB file has been successfully downloaded and saved to
            # the Logfiles directory
            if not os.path.exists(pathToPDB):
                sys.exit('Error 03: Failed to download and save PDB file - cause unknown')

        # If filepath to PDB has been supplied:
        else:
            # Changes directory to allow input PDB file (a '.pdb' or '.txt' file)
            # to be read from any provided file path. If PDB directory of this
            # name already exists in Logfiles directory, user input requested
            # ('Do you want to overwrite the existing file?'):
            # yes = old PDB directory is deleted, new PDB directory is created and
            #       copy of input PDB file is saved to the new directory
            # no = old PDB directory is retained, exit program
            # If it doesn't already exist, new PDB directory is created and copy of
            # input PDB file is saved to the new directory.
            owd = os.getcwd()
            pathToPDB = pathToPDB.replace('\\', '/')
            splitPath = pathToPDB.split('/')
            disk = '%s/' % splitPath[0]
            os.chdir('/')
            os.chdir(disk)

            if not os.path.exists(pathToPDB):
                sys.exit('Error 02: Supplied filepath not recognised')

            if pathToPDB[-4:] not in ['.pdb', '.txt']:
                sys.exit('Error 01: Supplied filepath to PDB is not a .pdb or'
                         '.txt file')
            else:
                print 'Filepath to .pdb or .txt file supplied\n'
                fileName = splitPath[len(splitPath)-1]
                splitFilename = fileName.split('.')
                fileName = splitFilename[len(splitFilename)-2].upper() + '.' + splitFilename[len(splitFilename)-1]
                PDBcode = splitFilename[len(splitFilename)-2].upper()
                PDBdirectory = 'Logfiles/%s/' % PDBcode
                pdb_file_path = '%s%s' % (PDBdirectory, PDBcode)
                newPathToPDB = '%s%s' % (PDBdirectory, fileName)

                os.chdir(owd)
                if os.path.isdir(PDBdirectory):
                    print 'Folder %s already exists locally at %s' % (PDBcode,
                                                                      PDBdirectory)
                    print('Do you want to overwrite the existing folder?\n'
                          '--USER INPUT-- type your choice and press RETURN\n'
                          'yes = overwrite this folder\n'
                          'no = do not overwrite this folder\n')
                    owChoice = None
                    while owChoice not in ['YES', 'NO', 'Y', 'N']:
                        owChoice = raw_input(prompt).upper()
                        if owChoice == 'YES' or owChoice == 'Y':
                            print '\nOverwriting existing folder'
                            shutil.rmtree(str(PDBdirectory))
                            copyPDB(pathToPDB, disk, newPathToPDB, PDBdirectory)
                            break
                        elif owChoice == 'NO' or owChoice == 'N':
                            print('\nKeeping original folder\n'
                                  'Exiting RABDAM')
                            sys.exit()
                            break
                        else:
                            print 'Unrecognised input - please answer "yes" or "no"'
                else:
                    copyPDB(pathToPDB, disk, newPathToPDB, PDBdirectory)

                # Checks that PDB file has been successfully copied to the Logfiles
                # directory
                if not os.path.exists(newPathToPDB):
                    sys.exit('Error 04: Failed to copy PDB to Logfiles directory.\n'
                             'Check that supplied PDB file is not in use by another'
                             'program')

                pathToPDB = newPathToPDB

        print('All files generated by this program will be stored in:\n'
              '    %s\n' % PDBdirectory)

        # Processes the input PDB file to remove hydrogen atoms, anisotropic
        # B factor records, and atoms with zero occupancy, as well as
        # retaining only the most probable alternate conformations. The unit cell
        # assembly is then generated from the coordinates of the processed PDB file
        # and its associated symmetry operations.
        print ('\nProcessing PDB file to remove hydrogen atoms, anisotropic '
               '\nB factor records, and atoms with zero occupancy, as well as '
               '\nretaining only the most probable alternate conformations')

        (clean_au_file, clean_au_list, header_lines, footer_lines,
         unit_cell_params) = clean_pdb_file(pathToPDB, pdb_file_path)

        print '\nGenerating unit cell\n'
        PDBCURinputFile = '%sPDBCURinput.txt' % pdb_file_path
        PDBCURlog = '%sPDBCURlog.txt' % pdb_file_path
        genPDBCURinputs(PDBCURinputFile)
        unit_cell_pdb = '%s_unit_cell.pdb' % pdb_file_path
        runPDBCUR(clean_au_file, unit_cell_pdb, PDBCURinputFile, PDBCURlog)

        if not os.path.exists(unit_cell_pdb):
            sys.exit('Error 05: Error in running PDBCUR, failed to generate Unit'
                     'Cell PDB file')

        print('****************** End of Process PDB Section ******************\n'
              '****************************************************************\n')

        print('****************************************************************\n'
              '********************** Parsing PDB Section *********************\n')
        # Parses the newly generated unit cell PDB file into RABDAM, returning a
        # list of all atoms in the unit cell, plus the unit cell parameters. The
        # unit cell PDB file is then deleted unless createUCpdb is set equal to
        # True in INPUT.txt.
        ucAtomList, xyzList = full_atom_list(unit_cell_pdb)

        if self.createUCpdb is False:
            os.remove(unit_cell_pdb)

        print('\n****************** End of Parsing PDB Section ******************\n'
              '****************************************************************\n')

        print('****************************************************************\n'
              '****************** Translate Unit Cell Section *****************\n')
        # The unit cell parameters are converted into Cartesian vectors. These
        # vectors are then used to translate the unit cell -/+ 1 units in all
        # (3) dimensions to generate a 3x3 parallelepiped. A PDB file of this 3x3
        # parallelepiped is output if createAUCpdb is set equal to True in
        # INPUT.txt.

        cartesianVectors = convertToCartesian(unit_cell_params)
        transAtomList = duplicate(ucAtomList)
        for a in xrange(-1, 2):
            for b in xrange(-1, 2):
                for c in xrange(-1, 2):
                    if a == 0 and b == 0 and c == 0:
                        pass  # No identity translation
                    else:
                        transAtomList = translateUnitCell(ucAtomList, transAtomList,
                                                          cartesianVectors, a, b, c)

        if self.createAUCpdb is True:
            aucPDBfilepath = '%s_all_unit_cells.pdb' % pdb_file_path
            makePDB(header_lines, transAtomList, footer_lines, aucPDBfilepath)

        print('\n************** End of Translate Unit Cell Section **************\n'
              '****************************************************************\n')

        print('****************************************************************\n'
              '********************* Trim Crystal Section *********************\n')
        # Atoms in the 3x3 parallelepiped which are a distance greater than the
        # packing density threshold from the boundaries of the processed
        # asymmetric unit are discarded. The PDB files of the processed
        # asymmetric unit and the trimmed parallelepiped are subsequently
        # deleted unless createAUpdb and createTApdb are respectively set equal
        # to True in INPUT.txt.
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

        print 'Removing atoms outside of packing density threshold\n'
        keepParams = convertParams(auParams, self.PDT)
        trimmedAtomList = trimAtoms(transAtomList, keepParams)

        if self.createTApdb is True:
            taPDBfilepath = '%s_trimmed_atoms.pdb' % pdb_file_path
            makePDB(header_lines, trimmedAtomList, footer_lines, taPDBfilepath)

        print('****************** End of Trim Crystal Section *****************\n'
              '****************************************************************\n')

        print('****************************************************************\n'
              '*************** Calculate Packing Density Section **************\n')
        # Calculates the packing density (atomic contact number) of all atoms
        # in the asymmetric unit.

        print 'Calculating packing density values\n'
        au_atom_xyz, trim_atom_xyz = get_xyz_from_objects(
            bdamAtomList, trimmedAtomList
            )
        packing_density_array = calc_packing_density(
            au_atom_xyz, trim_atom_xyz, self.PDT
            )
        write_pckg_dens_to_atoms(bdamAtomList, packing_density_array)

        print('*********** End of Calculate Packing Density Section ***********\n'
              '****************************************************************\n')

        print('****************************************************************\n'
              '****************** Calculate B_Damage Section ******************\n')
        # Atoms in the asymmetric unit are ordered via their packing density
        # values; the B_Damage value of each atom is then calculated as the ratio
        # of its B factor as compared to the average of the B factor values of
        # similarly (identified via sliding window) packed atoms.

        print 'Calculating B_Damage values\n'
        window = int(round((len(bdamAtomList)*self.windowSize), 0))
        if window % 2 == 0:
            window = window + 1  # Window size must be an odd number.
        if window < 3:
            window == 3  # Minimum window size is 3.
        calcBdam(bdamAtomList, window)

        print('****************************************************************\n'
              '******************* Writing DataFrame Section ******************\n')
        # Writes asymmetric unit atom properties, including their newly calculated
        # B_Damage values, to DataFrame, for ease of statistical analysis.
        # DataFrame, plus additional variables and lists required for subsequent
        # analysis, are pickled.

        print 'Writing B_Damage data to DataFrame\n'
        df = writeBdam(bdamAtomList)

        print 'Saving DataFrame\n'
        storage = '%s/DataFrame' % PDBdirectory
        os.makedirs(storage)
        storage_fileName = '%s/%s' % (storage, PDBcode)
        df.to_pickle(str(storage_fileName) + '_dataframe.pkl')
        with open(str(storage_fileName) + '_variables.pkl', 'wb') as f:
            pickle.dump((pdb_file_path, PDBcode, bdamAtomList, header_lines,
                         footer_lines, window), f)

        print('****************************************************************\n'
              '*************** End Of Writing DataFrame Section ***************\n')

        os.chdir('%s' % cwd)

    def rabdam_analysis(self, run):
        # Uses values in DataFrame returned from 'rabdam_dataframe' function to
        # write output analysis files.

        import os
        import sys
        prompt = '> '
        import pickle
        import pandas as pd
        from output import generate_output_files

        if run == 'rabdam_analysis':
            print '************************ RABDAM ANALYSIS ***********************\n'
            print('\nPlease cite: M. Gerstel, C. M. Deane and E.F. Garman. (2015).\n'
                  'J. Synchrotron Radiation, 22, 201-212\n'
                  'http://dx.doi.org/doi:10.1107/S1600577515002131\n')

        cwd = os.getcwd()
        os.chdir('%s' % self.outputDir)

        print('\n****************************************************************\n'
              '***************** Processing DataFrame Section *****************\n')
        # Checks that pkl files output by 'rabdam_dataframe' function required for
        # analysis exist.
        #
        # Then checks if Logfiles/PDBcode directory already contains any analysis
        # files - if so then user input requested:
        # yes = all old analysis files are removed and replaced by new analysis
        #       files
        # no = old analysis files are retained, exit program
        #
        # Note there is no option to replace only a subset of the output analysis
        # files.

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
        storage_fileName = '%s/%s' % (storage_directory, PDBcode)

        if not os.path.isdir(storage_directory):
            print 'Folder %s not found' % (storage_directory)
            print 'Exiting RABDAM analysis'
            sys.exit()

        potential_analysis_files = ['Bdamage.csv', 'Bdamage.html',
                                    'Bdamage.pdb', 'Bdamage.png',
                                    'Bnet_Protein.png', 'Bnet_NA.png']
        actual_analysis_files = []
        for name in potential_analysis_files:
            if os.path.isfile(str(PDB_analysis_file) + str(name)):
                actual_analysis_files.append(str(PDB_analysis_file) + str(name))
        if len(actual_analysis_files) > 0:
            print('There are one or more RABDAM analysis files already present\n'
                  'in folder %s' % PDBdirectory)
            print('Do you want to overwrite the existing analysis files?\n'
                  '--USER INPUT-- type your choice and press RETURN\n'
                  'yes = overwrite ALL analysis files\n'
                  'no = do not overwrite analysis files\n')
            owChoice = None
            while owChoice not in ['YES', 'NO', 'Y', 'N']:
                owChoice = raw_input(prompt).upper()
                if owChoice == 'YES' or owChoice == 'Y':
                    print 'Overwriting existing analysis files'
                    for name in actual_analysis_files:
                        os.remove(name)
                    break
                elif owChoice == 'NO' or owChoice == 'N':
                    print('Keeping original analysis files\n'
                          'Exiting RABDAM')
                    sys.exit()
                    break
                else:
                    print 'Unrecognised input - please answer "yes" or "no"'

        # Pkl files unpickled
        with open(str(storage_fileName) + '_variables.pkl', 'rb') as f:
            (pdb_file_path, PDBcode, bdamAtomList, header_lines, footer_lines,
             window) = pickle.load(f)
        df = pd.read_pickle(str(storage_fileName) + '_dataframe.pkl')

        print('************** End Of Processing DataFrame Section *************\n'
              '****************************************************************\n')

        print('****************************************************************\n'
              '***************** Writing Output Files Section *****************\n')
        # Uses values in DataFrame created by 'rabdam_dataframe' function to:
        # - draw a kernel density plot of distribution of B_Damage values of
        #   all analysed atoms in asymmetric unit
        # - write an html file listing all atoms whose B_Damage values lie above
        #   threshold specified in INPUT.txt
        # - calculate global B_Damage metric, plus draw a kernel density plot
        #   of B_Damage values of subset of atoms (Cys S, Glu O and Asp O)
        #   used to calculate it
        # - write pdb files with B factor values replaced by B_Damage values to
        #   allow structure when viewed with molecular graphics software to be
        #   coloured by B_Damage
        # - write DataFrame values to csv file, to allow user to easily further
        #   manipulate data as required
        output = generate_output_files(pdb_file_path=pdb_file_path, df=df)

        print 'Writing histogram and html files\n'
        output.make_histogram(self.threshold, self.highlightAtoms)

        print 'Calculating global B_Damage\n'
        output.calculate_Bnet(window_name, pdt_name)

        print 'Writing pdb files\n'
        output.make_colourbyBdam_pdb(header_lines, footer_lines, bdamAtomList)

        print 'Writing csv file\n'
        output.make_csv(bdamAtomList, window)

        print('************** End of Writing Output Files Section *************\n'
              '****************************************************************\n')

        os.chdir('%s' % cwd)
