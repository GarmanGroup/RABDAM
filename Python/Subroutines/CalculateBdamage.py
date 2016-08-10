

# Script to calculate B-damage for protein atoms
def rabdam(pathToPDB, PDT=14, binSize=10, HETATM=False, addAtoms=[], removeAtoms=[], threshold=0.02, createAUCpdb=False, createTApdb=False):
    print('\nRABDAM\n')
    print('\n')
    print('Please cite: M. Gerstel, C. M. Deane and E.F. Garman. (2015).\nJ. Synchrotron Radiation. 22, 201-212\nhttp://dx.doi.org/doi:10.1107/S1600577515002131\n')
    # import packages required for running the program
    import time  # for recording how long the script takes to run
    import sys  # for terminating script when encountering errors
    prompt = '> '
    import os  # for operating system usability
    import shutil
    import math  # for using more intricate mathematics
    import copy  # for making shallow copies of variables/lists/objects etc.
    duplicate = copy.copy
    from PDBCUR import genPDBCURinputs, runPDBCUR  # facilitates PDBCUR functionality
    from parsePDB import parsePDB, downloadPDB, copyPDB, getUnitCellParams, getAUparams, trimAtoms  # for taking information from PDB file to a usable format
    from translateUnitCell import convertToCartesian, getXYZlist, translateUnitCell  # translates unit cell
    from makePDB import makePDB, writeBdam  # allows new output files to be written from a list of atom objects
    from atomCheck import convertParams  # adds the PDT to the Cartesian limits of the unit cell
    from Bdamage import binAtoms, calcBdam, get_xyz_from_objects, calc_packing_density, write_pckg_dens_to_atoms  # calculates the PDT of each atom in the proivided structure
    from output import make_csv, make_histogram, make_colourbyBdam_pdb
    # Input: the file path to the pdb for which you want to calculate B-damage factors, the 'Packing Density Threshold' (Angstroms) and bin size
    start = time.time()
    startIndex = time.gmtime()
    year = startIndex.tm_year
    month = startIndex.tm_mon
    day = startIndex.tm_mday
    hour = startIndex.tm_hour
    minute = startIndex.tm_min
    second = startIndex.tm_sec
    # inform the user of the start time of the program
    print 'This program was run on %d/%d/%d at %02.0f:%02.0f:%02.0f\n\n' % (day, month, year, hour, minute, second)
    print '################################################################'
    print '################################################################'
    print '########## Program to calculate B Damage #######################'
    print '################################################################'
    print '\n'

    # inform the user of the inputs used in the program
    print '****************************************************************'
    print '********** Input Section ***************************************\n'
    print 'Calculating B Damage for %s' % pathToPDB
    if PDT == 14:
        print 'No packing density threshold supplied, using default value of\n14 Angstroms\n'
    else:
        print 'Packing density threshold defined by user at %s Angstroms\n' % PDT
    if binSize == 10:
        print 'No bin size supplied, using default value of 10\n'
    else:
        print 'Bin size defined by user as %s\n' % binSize
    if HETATM is True:
        print 'Waters to be removed, other HETATM to be retained\n'
    elif HETATM is False:
        print 'All HETATM to be removed\n'
    if len(addAtoms) == 0:
        print 'No atoms to be added\n'
    else:
        for index, value in enumerate(addAtoms):
            print 'Atoms to be added: %s\n' % addAtoms[index]
    if len(removeAtoms) == 0:
        print 'No atoms to be removed\n'
    else:
        for index, value in enumerate(removeAtoms):
            print 'Atoms to be removed: %s\n' % removeAtoms[index]
    if threshold == 0.02:
        print 'No threshold value supplied, using default value of 0.02\n'
    else:
        print 'Threshold value defined by user as %s\n' % threshold
    print '********** End of Input Section ********************************'
    print '****************************************************************'
    print '\n'
    # Process pdb file with PDBCUR to generate a unit cell and then store the name of the processed pdb
    print '****************************************************************'
    print '********** Process PDB Section *********************************\n'
    print 'Processing PDB file to:'
    print '- Remove Hydrogen atoms'
    print '- Remove atoms with zero occupancy.'
    print '- Keep most probable alternate conformation.'
    print '- Remove anisotropic B factor records from the file.'
    print '- Generate the Unit Cell from Symmetry Operations.\n'
    # check whether filepath or PDB code supplied
    if len(pathToPDB) == 4:
        print 'PDB code supplied'
        # convert PDB accession code to UPPERCASE
        PDBcode = pathToPDB.upper()
        PDBdirectory = 'Logfiles/%s/' % PDBcode
        pathToPDB = '%s%s.pdb' % (PDBdirectory, PDBcode)
        # Check if PDB file has already been downloaded
        if os.path.isdir(PDBdirectory):
            # inform user that file already exists
            print 'Folder %s already exists locally at %s' % (PDBcode, PDBdirectory)
            print 'Do you want to overwrite the existing file?\n'
            print '--USER INPUT-- type your choice and press RETURN\n'
            print 'yes = overwrite this folder'
            print 'no = do not overwrite this folder'
            owChoice = None
            while owChoice not in ['yes', 'no', 'Yes', 'No', 'YES', 'NO', 'y', 'n', 'Y', 'N']:
                owChoice = raw_input(prompt)
                if owChoice == 'yes' or owChoice == 'Yes' or owChoice == 'YES' or owChoice == 'y' or owChoice == 'Y':
                    print 'overwriting existing folder'
                    shutil.rmtree(str(PDBdirectory))
                    downloadPDB(PDBcode, PDBdirectory, pathToPDB)
                    break
                elif owChoice == 'no' or owChoice == 'No' or owChoice == 'NO' or owChoice == 'n' or owChoice == 'N':
                    print 'keeping original folder'
                    print 'exiting RABDAM'
                    sys.exit()
                    break
                else:
                    print 'unrecognised input - please answer yes or no'
        else:
            downloadPDB(PDBcode, PDBdirectory, pathToPDB)
            owChoice = 'null'
            # check that file has downloaded and saved correctly
        if not os.path.exists(pathToPDB):
            sys.exit('Error 03: Failed to download and save PDB - cause unknown')
    else:
        # change directory to allow file to be read from any path on local disk
        owd = os.getcwd()
        os.chdir('c://')
        # check supplied filepath exists, returning error message if not
        if not os.path.exists(pathToPDB):
            sys.exit('Error 02: Supplied filepath does not exist')
        # check supplied filepath is a pdb file, returning error message if not
        if pathToPDB[-4:] == '.pdb' or pathToPDB[-4:] == '.txt':
            print 'Filepath to .pdb or .txt file supplied'
            # copy file to Logfiles directory if necessary
            splitPath = pathToPDB.split('/')
            fileName = splitPath[len(splitPath)-1]
            splitFilename = fileName.split(".")
            PDBcode = splitFilename[len(splitFilename)-2]
            PDBdirectory = 'Logfiles/%s/' % PDBcode
            newPathToPDB = '%s%s' % (PDBdirectory, fileName)
            # check if file has already been copied to Logfiles
            os.chdir(owd)
            if os.path.isdir(PDBdirectory):
                # inform the user that file already exists
                print 'Folder %s already exists locally at %s' % (PDBcode, PDBdirectory)
                print 'Do you want to overwrite the existing folder?\n'
                print '--USER INPUT-- type your choice and press RETURN\n'
                print 'yes = overwrite this folder'
                print 'no = do not overwrite this folder'
                owChoice = None
                while owChoice not in ['yes', 'no', 'Yes', 'No', 'YES', 'NO', 'y', 'n', 'Y', 'N']:
                    owChoice = raw_input(prompt)
                    if owChoice == 'yes' or owChoice == 'Yes' or owChoice == 'YES' or owChoice == 'y' or owChoice == 'Y':
                        print 'overwriting existing folder'
                        shutil.rmtree(str(PDBdirectory))
                        copyPDB(pathToPDB, newPathToPDB, PDBdirectory)
                        break
                    elif owChoice == 'no' or owChoice == 'No' or owChoice == 'NO' or owChoice == 'n' or owChoice == 'N':
                        print 'keeping original folder'
                        print 'exiting RABDAM'
                        sys.exit()
                        break
                    else:
                        print 'unrecognised input - please answer yes or no'
            else:
                # make local copy in Logfiles
                copyPDB(pathToPDB, newPathToPDB, PDBdirectory)
                owChoice = 'null'
            # check that file has copied and saved correctly
            if not os.path.exists(newPathToPDB):
                sys.exit('Error 04: Failed to copy PDB to a local version.\nCheck that supplied PDB is not in use by another program')
            # rename pathToPDB variable so that files generated by RABDAM are saved in the PDB code folder
            pathToPDB = newPathToPDB
        else:
            sys.exit('Error 01: Supplied filepath to PDB is not a .pdb or .txt file')

    print 'All files generated by this program will be stored in \n%s\n' % PDBdirectory
    # create path to PDBCUR input file
    PDBCURinputFile = '%sPDBCURinput.txt' % PDBdirectory
    # generate path to PDBCUR log file
    PDBCURlog = '%sPDBCURlog.txt' % PDBdirectory
    # generate input file for PDBCUR
    genPDBCURinputs(PDBCURinputFile)
    # Create name for output PDBCUR file by appending 'UnitCell' to filename
    splitFilePath = pathToPDB.split(".")
    fileName = splitFilePath[len(splitFilePath)-2]
    PDBCURoutputPDB = '%sUnitCell.pdb' % fileName
    # runPDBCUR using generated input file
    runPDBCUR(pathToPDB, PDBCURoutputPDB, PDBCURinputFile, PDBCURlog)
    if not os.path.exists(PDBCURoutputPDB):
        sys.exit('Error 05: Failed to generate Unit Cell PDB file')
    print '\n********** End of Process PDB Section **************************'
    print '****************************************************************'
    print '\n'
    # Parse the newly generated unit cell PDB file
    print '****************************************************************'
    print '********** Parsing PDB Section *********************************\n'
    # return a list of atoms and attributes
    bof, ucAtomList, bdamAtomList, eof = parsePDB(PDBCURoutputPDB, HETATM, addAtoms, removeAtoms)
    unitCell = getUnitCellParams(pathToPDB)
    print '\n********** End of Parsing PDB Section **************************'
    print '****************************************************************'
    print '\n'
    # Translate the Unit Cell and append the new atom locations to the ucAtomList
    print '****************************************************************'
    print '********** Translate Unit Cell Section *************************\n'
    # convert the unit cell paramenters to Cartesian coordinates
    cartesianVectors = convertToCartesian(unitCell)
    # obtain an array of XYZcoordinates from input list of atom objects
    xyzList = getXYZlist(ucAtomList)
    # create shallow copy of the list of atoms to which all translated atomic positions will be added
    transAtomList = duplicate(ucAtomList)
    # loop through running the translation subroutine for all combinations of
    # translations +/- 1 unit cell in a, b and c directions
    for a in xrange(-1, 2):
        for b in xrange(-1, 2):
            for c in xrange(-1, 2):
                # don't translate cells ine the case of an identity translation
                if a == 0 and b == 0 and c == 0:
                    pass
                else:
                    transAtomList = translateUnitCell(xyzList, ucAtomList, transAtomList, cartesianVectors, a, b, c)

    print ''
    if createAUCpdb is True:
        aucPDBfilepath = '%sAllUnitCells.pdb' % fileName
        makePDB(bof, transAtomList, eof, aucPDBfilepath)
    print '\n********** End of Translate Unit Cell Section ******************'
    print '****************************************************************'
    print '\n'
    # Discard atoms too far from the asymmetric unit
    print '****************************************************************'
    print '********** Trim Crystal Section ********************************\n'
    # parse a new set of atomic coordinates from the provided asymmetric unit file
    bof1, auAtomList, bdamAtomList, eof1 = parsePDB(pathToPDB, HETATM, addAtoms, removeAtoms)
    bof1.remove
    eof1.remove
    auParams = getAUparams(auAtomList)
    print 'Obtained asymmetric unit parameters:'
    print 'xMin = %8.3f' % auParams[0]
    print 'xMax = %8.3f' % auParams[1]
    print 'yMin = %8.3f' % auParams[2]
    print 'yMax = %8.3f' % auParams[3]
    print 'zMin = %8.3f' % auParams[4]
    print 'zMax = %8.3f\n' % auParams[5]
    # add PDT to auParams in all dimensions
    print 'Now creating a box surrounding the atoms'
    print 'We only want to consider atoms in this box when we calculate the'
    print 'packing density'
    keepParams = convertParams(auParams, PDT)
    # discard all atoms not in cube defined above
    trimmedAtomList = trimAtoms(transAtomList, keepParams)
    # create PDB file of retained atoms
    if createTApdb is True:
        taPDBfilepath = '%sTrimmedAtoms.pdb' % fileName
        makePDB(bof, trimmedAtomList, eof, taPDBfilepath)
    print '\n********** End of Trim Crystal Section *************************'
    print '****************************************************************'
    print '\n'
    # Calculate the packing density of each atom
    print '****************************************************************'
    print '********** Calculate Packing Density Section *******************\n'
    au_atom_xyz, trim_atom_xyz = get_xyz_from_objects(bdamAtomList, trimmedAtomList)
    packing_density_array = calc_packing_density(au_atom_xyz, trim_atom_xyz, PDT)
    minPD, maxPD = int(packing_density_array.min()), int(packing_density_array.max())
    write_pckg_dens_to_atoms(bdamAtomList, packing_density_array)
    noOfGroups, adjtNo = binAtoms(bdamAtomList, binSize, minPD, maxPD)
    print 'Lowest calculated Packing Density was %s' % minPD
    print 'Highest calculated Packing Density was %s' % maxPD
    print 'Atoms separated into %s bins\n' % noOfGroups
    print '\n********** End of Calculate Packing Density Section ************'
    print '****************************************************************'
    print '\n'
    # Calculate the Bdamage value of each atom
    print '****************************************************************'
    print '********** Calculate Bdamage Section ***************************\n'
    print 'Calculating Bdamage values'
    groupNoAtoms, groupAvBfacs = calcBdam(bdamAtomList, noOfGroups)
    print 'Writing Bdamage data to output file'

    df = writeBdam(groupNoAtoms, groupAvBfacs, binSize, minPD, adjtNo, bdamAtomList)
    x_values_RHS = make_histogram(df, fileName, PDBcode, threshold)

    bDamFileName = '%sBdamage.csv' % fileName
    make_csv(bdamAtomList, bDamFileName, groupNoAtoms, groupAvBfacs, binSize, adjtNo)
    make_colourbyBdam_pdb(df, bof, eof, fileName, bdamAtomList, x_values_RHS)

    print '\n********** End of Calculate Bdamage Section ********************'
    print '****************************************************************'
    print '\n'

    # inform the user of the time elapsed while the program was run
    runtime = time.time() - start
    minutes = math.floor(runtime/60)
    seconds = math.fmod(runtime, 60)
    if minutes == 0:
        if seconds == 1:
            print 'Total time taken for program to run was %02.3f second.\n\n' % seconds
        else:
            print 'Total time taken for program to run was %02.3f seconds.\n\n' % seconds
    elif minutes == 1:
        print 'Total time taken for program to run was %01.0f minute and %02.3f seconds.\n\n' % (minutes, seconds)
    else:
        print 'Total time taken for program to run was %01.0f minutes and %02.3f seconds.\n\n' % (minutes, seconds)
# end
