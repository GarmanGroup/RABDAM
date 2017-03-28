

def rabdam_dataframe(pathToPDB, PDT=14, windowSize=0.02,
                     protOrNA='PROTEIN', HETATM=False, addAtoms=[],
                     removeAtoms=[], createAUpdb=False, createUCpdb=False,
                     createAUCpdb=False, createTApdb=False, run='rabdam'):
    # Calculates B_damage values from input PDB file and writes
    # output to DataFrame.

    import sys
    prompt = '> '
    import os
    import shutil
    import copy
    duplicate = copy.copy
    import pickle

    from PDBCUR import genPDBCURinputs, runPDBCUR
    from parsePDB import (full_atom_list, b_damage_atom_list, downloadPDB,
                          copyPDB, getUnitCellParams, getAUparams, trimAtoms)
    from translateUnitCell import (convertToCartesian, getXYZlist,
                                   translateUnitCell)
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
          '**************** Program to calculate B_damage *****************\n'
          '****************************************************************\n')

    print('****************************************************************\n'
          '************************* Input Section ************************\n')
    # Prints the input values read into the program from INPUT.txt.

    print 'Calculating B_damage for %s\n' % pathToPDB
    if PDT == 14:
        print 'Using default packing density threshold of 14 Angstroms\n'
    else:
        print 'Packing density threshold defined by user as %s Angstroms\n' % PDT
    if windowSize == 0.02:
        print 'Using default window size of 0.02\n'
    else:
        print 'Window size defined by user as %s\n' % windowSize
    if HETATM is True:
        print 'All HETATM to be retained\n'
    elif HETATM is False:
        print 'All HETATM to be removed\n'
    if protOrNA == 'PROTEIN':
        print 'Retaining protein atoms, discarding nucleic acid atoms\n'
    elif protOrNA in ['NUCLEICACID', 'NA']:
        print 'Retaining nucleic acid atoms, discarding protein atoms\n'
    if len(addAtoms) == 0:
        print 'No atoms to be added\n'
    else:
        for value in addAtoms:
            print 'Atoms to be added: %s\n' % value
    if len(removeAtoms) == 0:
        print 'No atoms to be removed\n'
    else:
        for value in removeAtoms:
            print 'Atoms to be removed: %s\n' % value

    print('********************* End of Input Section *********************\n'
          '****************************************************************\n')

    print('****************************************************************\n'
          '********************** Process PDB Section *********************\n')
    # Creates a new directory named after the input PDB file (after checking
    # that this directory does not already exist) in the Logfiles directory.
    # Then saves a copy of the input PDB file to the new directory.

    # If 4 digit PDB accession code has been supplied:
    if len(pathToPDB) == 4:
        print 'PDB code supplied\n'
        PDBcode = pathToPDB.upper()
        PDBdirectory = 'Logfiles/%s/' % PDBcode
        pathToPDB = '%s%s.pdb' % (PDBdirectory, PDBcode)

        # If PDB directory of this name already exists in Logfiles directory,
        # user input requested ('Do you want to overwrite the existing file?'):
        # yes = old PDB directory is deleted, new PDB directory is created and
        #       copy of the PDB file is downloaded from the RSCB PDB website
        #       and saved to the new directory
        # no = old PDB directory is retained, exit program
        # If doesn't already exist, new PDB directory is created and copy of
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
            owChoice = 'null'

        if not os.path.exists(pathToPDB):
            sys.exit('Error 03: Failed to download and save PDB - cause unknown')

    # If filepath to PDB has been supplied:
    else:
        # Changes directory to allow input PDB file (a '.pdb' or '.txt' file)
        # to be read from any path on local disk (C://). If PDB directory of
        # this name already exists in Logfiles directory, user input
        # requested ('Do you want to overwrite the existing file?'):
        # yes = old PDB directory is deleted, new PDB directory is created and
        #       copy of input PDB file is saved to the new directory
        # no = old PDB directory is retained, exit program
        # If doesn't already exist, new PDB directory is created and copy of
        # input PDB file is saved to the new directory.
        owd = os.getcwd()
        pathToPDB = pathToPDB.replace('\\', '/')
        splitPath = pathToPDB.split('/')
        disk = '%s/' % splitPath[0]
        os.chdir(disk)

        if not os.path.exists(pathToPDB):
            sys.exit('Error 02: Supplied filepath not recognised')

        if pathToPDB[-4:] == '.pdb' or pathToPDB[-4:] == '.txt':
            print 'Filepath to .pdb or .txt file supplied\n'
            fileName = splitPath[len(splitPath)-1]
            splitFilename = fileName.split('.')
            fileName = str((splitFilename[len(splitFilename)-2]).upper()) + '.' + str(splitFilename[len(splitFilename)-1])
            PDBcode = splitFilename[len(splitFilename)-2].upper()
            PDBdirectory = 'Logfiles/%s/' % PDBcode
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
                owChoice = 'null'

            if not os.path.exists(newPathToPDB):
                sys.exit('Error 04: Failed to copy PDB to Logfiles directory.\n'
                         'Check that supplied PDB is not in use by another'
                         'program')

            pathToPDB = newPathToPDB

        else:
            sys.exit('Error 01: Supplied filepath to PDB is not a .pdb or'
                     '.txt file')

    print('All files generated by this program will be stored in:'
          '\n%s\n' % PDBdirectory)

    # Processes the input PDB file via PDBCUR to remove hydrogen atoms,
    # atoms with zero occupancy and anisotropic B factor records, and to
    # retain only the most probable alternate conformations and set their
    # occupancies to 1. PDBCUR is then used to generate a unit cell from the
    # processed PDB file and its associated symmetry operations.

    print '\nProcessing PDB file with PDBCUR\n'

    splitFilePath = pathToPDB.split('.')
    fileName = splitFilePath[len(splitFilePath)-2]
    asymmetricUnit = False  # Instructs PDBCUR to generate unit cell from asymmetric unit.
    PDBCURinputFile = '%sPDBCURinput.txt' % fileName
    PDBCURlog = '%sPDBCURlog.txt' % fileName
    genPDBCURinputs(PDBCURinputFile, asymmetricUnit)
    PDBCURoutputPDBuc = '%sUnitCell.pdb' % fileName
    runPDBCUR(pathToPDB, PDBCURoutputPDBuc, PDBCURinputFile, PDBCURlog)

    if not os.path.exists(PDBCURoutputPDBuc):
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

    bof, ucAtomList, eof = full_atom_list(PDBCURoutputPDBuc)

    if createUCpdb is False:
        os.remove(PDBCURoutputPDBuc)

    unitCell = getUnitCellParams(pathToPDB)

    print('\n****************** End of Parsing PDB Section ******************\n'
          '****************************************************************\n')

    print('****************************************************************\n'
          '****************** Translate Unit Cell Section *****************\n')
    # The unit cell parameters are converted into Cartesian vectors. These
    # vectors are then used to translate the unit cell -/+ 1 units in all
    # (3) dimensions to generate a 3x3 parallelepiped. A PDB file of this 3x3
    # parallelepiped is output if createAUCpdb is set equal to True in
    # INPUT.txt.

    cartesianVectors = convertToCartesian(unitCell)
    xyzList = getXYZlist(ucAtomList)
    transAtomList = duplicate(ucAtomList)
    for a in xrange(-1, 2):
        for b in xrange(-1, 2):
            for c in xrange(-1, 2):
                if a == 0 and b == 0 and c == 0:
                    pass  # No identity translation
                else:
                    transAtomList = translateUnitCell(
                        xyzList, ucAtomList, transAtomList, cartesianVectors,
                        a, b, c
                        )

    if createAUCpdb is True:
        aucPDBfilepath = '%sAllUnitCells.pdb' % fileName
        makePDB(bof, transAtomList, eof, aucPDBfilepath)

    print('\n************** End of Translate Unit Cell Section **************\n'
          '****************************************************************\n')

    print('****************************************************************\n'
          '********************* Trim Crystal Section *********************\n')
    # PDBCUR is run to process the asymmetric unit to remove hydrogen atoms,
    # atoms with zero occupancy and anisotropic B factor records, and to
    # retain only the most probable alternate conformations and set their
    # occupancies to 1. Atoms in the 3x3 parallelepiped which are a distance
    # greater than the packing density threshold from the boundaries of the
    # processed asymmetric unit are then discarded. The PDB files of the
    # processed asymmetric unit and the trimmed parallelepiped are subsequently
    # deleted unless createAUpdb and createTApdb are respectively set equal to
    # True in INPUT.txt.

    asymmetricUnit = True  # Instructs PDBCUR not to generate unit cell from asymmetric unit.
    PDBCURinputFile = '%sPDBCURinput.txt' % fileName
    PDBCURlog = '%sPDBCURlog.txt' % fileName
    genPDBCURinputs(PDBCURinputFile, asymmetricUnit)
    PDBCURoutputPDBau = '%sAsymmetricUnit.pdb' % fileName
    runPDBCUR(pathToPDB, PDBCURoutputPDBau, PDBCURinputFile, PDBCURlog)

    if not os.path.exists(PDBCURoutputPDBau):
        sys.exit('Error 05: Error in running PDBCUR, failed to generate'
                 'Asymmetric Unit PDB file')

    bof1, auAtomList, eof1 = full_atom_list(PDBCURoutputPDBau)
    bof1.remove
    eof1.remove
    bdamAtomList = b_damage_atom_list(
        PDBCURoutputPDBau, auAtomList, HETATM, protOrNA, addAtoms, removeAtoms
        )

    if createAUpdb is False:
        os.remove(PDBCURoutputPDBau)

    auParams = getAUparams(bdamAtomList)
    print '\nObtained asymmetric unit parameters:'
    print 'xMin = %8.3f' % auParams[0]
    print 'xMax = %8.3f' % auParams[1]
    print 'yMin = %8.3f' % auParams[2]
    print 'yMax = %8.3f' % auParams[3]
    print 'zMin = %8.3f' % auParams[4]
    print 'zMax = %8.3f\n' % auParams[5]

    print 'Removing atoms outside of packing density threshold\n'
    keepParams = convertParams(auParams, PDT)
    trimmedAtomList = trimAtoms(transAtomList, keepParams)

    if createTApdb is True:
        taPDBfilepath = '%sTrimmedAtoms.pdb' % fileName
        makePDB(bof, trimmedAtomList, eof, taPDBfilepath)

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
        au_atom_xyz, trim_atom_xyz, PDT
        )
    write_pckg_dens_to_atoms(bdamAtomList, packing_density_array)

    print('*********** End of Calculate Packing Density Section ***********\n'
          '****************************************************************\n')

    print('****************************************************************\n'
          '****************** Calculate B_damage Section ******************\n')
    # Atoms in the asymmetric unit are ordered via their packing density
    # values; the B_damage value of each atom is then calculated as the ratio
    # of its B factor as compared to the average of the B factor values of
    # similarly (identified via sliding window) packed atoms.

    print 'Calculating B_damage values\n'
    window = int(round((len(bdamAtomList)*windowSize), 0))
    if window % 2 == 0:
        window = window + 1  # Window size must be an odd number.
    if window < 15:
        window == 15  # Minimum window size is 15.
    calcBdam(bdamAtomList, window)

    print('****************************************************************\n'
          '******************* Writing DataFrame Section ******************\n')
    # Writes asymmetric unit atom properties, including their newly calculated
    # B_damage values, to DataFrame, for ease of statistical analysis.
    # DataFrame, plus additional variables and lists required for subsequent
    # analysis, are pickled.

    print 'Writing B_damage data to DataFrame\n'
    df = writeBdam(bdamAtomList)

    print 'Saving DataFrame\n'
    storage = '%s/DataFrame' % PDBdirectory
    os.makedirs(storage)
    storage_fileName = '%s/%s' % (storage, PDBcode)
    df.to_pickle(str(storage_fileName) + '_dataframe.pkl')
    with open(str(storage_fileName) + '_variables.pkl', 'wb') as f:
        pickle.dump((fileName, PDBcode, bdamAtomList, bof, eof, window), f)

    print('****************************************************************\n'
          '*************** End Of Writing DataFrame Section ***************\n')


def rabdam_analysis(pathToPDB, threshold=0.02, highlightAtoms=[],
                    run='rabdam'):
    # Uses values in DataFrame returned from 'rabdam_dataframe' function to
    # write output analysis files.

    import os
    import sys
    prompt = '> '
    import pickle
    import pandas as pd
    from output import (make_csv, make_histogram, make_colourbyBdam_pdb)

    if run == 'rabdam_analysis':
        print '************************ RABDAM ANALYSIS ***********************\n'
        print('\nPlease cite: M. Gerstel, C. M. Deane and E.F. Garman. (2015).\n'
              'J. Synchrotron Radiation, 22, 201-212\n'
              'http://dx.doi.org/doi:10.1107/S1600577515002131\n')

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

    splitPath = pathToPDB.split('/')
    pathToPDB = splitPath[len(splitPath)-1]
    PDBcode = pathToPDB.replace('.pdb', '')
    PDBcode = PDBcode.replace('.txt', '')
    PDBcode = PDBcode.upper()
    PDBdirectory = 'Logfiles/%s' % PDBcode
    PDB_analysis_file = '%s/%s' % (PDBdirectory, PDBcode)
    storage_directory = '%s/DataFrame' % PDBdirectory
    storage_fileName = '%s/%s' % (storage_directory, PDBcode)

    if not os.path.isdir(storage_directory):
        print 'Folder %s not found' % (storage_directory)
        print 'Exiting RABDAM analysis'
        sys.exit()

    potential_analysis_files = [
                 'Bdamage.csv', 'Bdamage.html', 'Bdamage.pdb',
                 'Bdamage.png', 'Bdamage_above_boundary.pdb'
                 ]
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
        fileName, PDBcode, bdamAtomList, bof, eof, window = pickle.load(f)
    df = pd.read_pickle(str(storage_fileName) + '_dataframe.pkl')

    print('************** End Of Processing DataFrame Section *************\n'
          '****************************************************************\n')

    print('****************************************************************\n'
          '***************** Writing Output Files Section *****************\n')
    # Uses values in DataFrame created by 'rabdam_dataframe' function to:
    # - draw a kernel density plot of distribution of B_damage values of
    #   all analysed atoms in asymmetric unit
    # - write an html file listing all atoms whose B_damage values lie above
    #   threshold specified in INPUT.txt
    # - calculate global B_damage metric, plus draw a kernel density plot
    #   of B_damage values of subset of atoms (Cys S, Glu O and Asp O)
    #   used to calculate it
    # - write pdb files with B factor values replaced by B_damage values to
    #   allow structure when viewed with molecular graphics software to be
    #   coloured by B_damage
    # - write DataFrame values to csv file, to allow user to easily further
    #   manipulate data as required

    print 'Writing histogram and html files\n'
    x_values_RHS = make_histogram(
        df, fileName, PDBcode, threshold, highlightAtoms
        )

    print 'Calculating global B_damage\n'
    #  calculate_global_BDam(df, PDBcode, fileName)

    print 'Writing pdb files\n'
    make_colourbyBdam_pdb(df, bof, eof, fileName, bdamAtomList, x_values_RHS)

    print 'Writing csv file\n'
    bDamFileName = '%sBdamage.csv' % fileName
    make_csv(bdamAtomList, bDamFileName, window)

    print('************** End of Writing Output Files Section *************\n'
          '****************************************************************\n')
