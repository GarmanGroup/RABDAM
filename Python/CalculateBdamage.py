# Copyright 2015 Thomas Dixon
# With thanks to Jonathan Brooks-Bartlett, Charles Bury, Markus Gerstel and Elspeth Garman

#Script to calculate B-damage for protein atoms
def CalculateBdamage(pathToPDB, PDT=14, binSize=10, createAllUnitCellsPDB=1, createTrimmedAtomsPDB=1):
    print('\n')
    print('Copyright 2015 Thomas Dixon\n')
    print('With thanks to Jonathan Brooks-Bartlett, Charles Bury, Markus\nGerstel and Elspeth Garman')
    #import packages required for running the program
    import time #for recording how long the script takes to run
    import sys #for terminating script when encountering errors
    import urllib2 #for dealing with URL stuff
    import os #for operating system usability
    import math #for using more intricate mathematics
    import copy #for making shallow copies of variables/lists/objects etc.
    from PDBCUR import genPDBCURinputs,runPDBCUR #facilitates PDBCUR functionality
    from parsePDB import parsePDB, getUnitCellParams, getAUparams, trimAtoms #for taking information from PDB file to a usable format
    from translateUnitCell import convertToCartesian,translateUnitCell #translates unit cell
    from makePDB import makePDB, writeBdam #allows new PDB files to be written from a list of atom objects
    from atomCheck import convertParams #adds the PDT to the Cartesian limits of the unit cell
    from Bdamage import calcPDT, binAtoms, calcBdam #calculates the PDT of each atom in the proivided structure
    #Input: the file path to the pdb for which you want to calculate B-damage factors, the 'Packing Density Threshold' (Angstroms) and bin size
    start = time.time()
    startIndex = time.gmtime()
    year = startIndex.tm_year
    month = startIndex.tm_mon
    day = startIndex.tm_mday
    hour = startIndex.tm_hour
    minute = startIndex.tm_min
    second = startIndex.tm_sec
    #inform the user of the start time of the program
    print 'This program was run on %d/%d/%d at %02.0f:%02.0f:%02.0f\n\n' % (day, month, year, hour, minute, second)
    print '################################################################'
    print '################################################################'
    print '########## Program to calculate B Damage #######################'
    print '################################################################'
    print '\n'
    #inform the user of the inputs used in the program
    print '****************************************************************'
    print '********** Input Section ***************************************\n'
    print 'Calculating B Damage for %s' % pathToPDB
    if PDT == 14:
        print 'No packing density threshold supplied, using default value of\n14 Angstroms'
    else:
        print 'Packing density threshold defined by user at %s Angstroms\n' % PDT
    if binSize == 10:
        print 'No bin size supplied, using default value of 10\n'
    else:
        print 'Bin size defined by user as %s\n' % binSize
    print '********** End of Input Section ********************************'
    print '****************************************************************'
    print '\n'
    #Process pdb file with PDBCUR to generate a unit cell and then store the name of the processed pdb
    print '****************************************************************'
    print '********** Process PDB Section *********************************\n'
    print 'Processing PDB file to:'
    print '- Remove Hydrogen atoms'
    print '- Remove atoms with zero occupancy.'
    print '- Keep most probable alternate conformation.'
    print '- Remove anisotropic B factor records from the file.'
    print '- Generate the Unit Cell from Symmetry Operations.\n'
    #check whether filepath or PDB code supplied
    if len(pathToPDB) == 4:
        print 'PDB code supplied'
        #convert PDB accession code to UPPERCASE
        PDBcode = pathToPDB.upper()
        PDBdirectory = 'Logfiles/%s/' % PDBcode
        pathToPDB = '%s%s.pdb' % (PDBdirectory, PDBcode)
        #Check if PDB file has already been downloaded
        if os.path.isfile(pathToPDB):
            #inform user that file already exists
            print 'PDB file already exists locally at %s' % pathToPDB
        else:
            #create URL from which to download .pdb file
            urlText = 'http://www.rcsb.org/pdb/files/%s.pdb' % PDBcode
            #downlaod PDB file and save local copy
            if os.path.exists(PDBdirectory):
                print 'Directory %s already exists' % PDBdirectory
            else:
                os.makedirs(PDBdirectory)
                print 'Directory %s created' % PDBdirectory
            origPDB = urllib2.urlopen(urlText)
            #inform user of the URL used to download PDB file
            print 'Downloaded PDB file from %s' % urlText
            #write local file containing the downloaded content
            localFile = open(pathToPDB, 'w') 
            localFile.write(origPDB.read())
            #inform user of file loaction of newly downloaded content
            print 'PDB file saved to %s' % pathToPDB
            #close local file to free up memory
            localFile.close()
            #chack that file has downloaded and saved correctly
            if not os.path.exists(pathToPDB):
                sys.exit ('Error 03: Failed to download and save PDB - cause unknown')
    else:
        #check supplied filepath is a pdb file, returning error message if not
        if pathToPDB[-4:] == '.pdb':
            print 'Filepath to .pdb file supplied'
            #copy file to Logfiles directory if necessary
            splitPath = pathToPDB.split("/")
            fileName = splitPath[len(splitPath)-1]
            splitFilename = fileName.split(".")
            directoryName = splitFilename[len(splitFilename)-2]
            PDBdirectory = 'Logfiles/%s/' % directoryName
            newPathToPDB = '%s%s' % (PDBdirectory, fileName)
            #check if file has already been copied to Logfiles
            if os.path.isfile(newPathToPDB):
                #inform the user that file already exists
                print 'PDB file already exists locally at %s' % newPathToPDB
            else:
                #make local copy in Logfiles
                if not os.path.exists(PDBdirectory):
                    os.makedirs(PDBdirectory)
                origPDB = open(pathToPDB, 'r')
                #write file containing copied content
                localFile = open(newPathToPDB, 'w')
                localFile.write(origPDB.read())
                #inform user of file location of new copy of PDB file
                print 'PDB file copied to %s' % newPathToPDB
                #close local file to free up memory
                localFile.close()
            #chack that file has downloaded and saved correctly
            if not os.path.exists(newPathToPDB):
                sys.exit ('Error 04: Failed to copy PDB to a local version. Check that supplied PDB is not in use by another program')
        else:
            sys.exit('Error 01: Supplied filepath to PDB is not a .pdb file')
        #check supplied filepath exists, returning error message if not
        if not os.path.exists(pathToPDB):
            sys.exit ('Error 02: Supplied filepath does not exist')
    print 'All files generated by this program will be stored in \n%s\n' % PDBdirectory
    #create path to PDBCUR input file
    PDBCURinputFile = '%sPDBCURinput.txt' % PDBdirectory
    #generate path to PDBCUR log file
    PDBCURlog = '%sPDBCURlog.txt' % PDBdirectory
    #generate input file for PDBCUR
    genPDBCURinputs(PDBCURinputFile)
    #Create name for output PDBCUR file by appending 'UnitCell' to filename 
    splitFilePath = pathToPDB.split(".")
    fileName = splitFilePath[len(splitFilePath)-2]
    PDBCURoutputPDB = '%sUnitCell.pdb' % fileName
    #runPDBCUR using generated input file
    runPDBCUR(pathToPDB, PDBCURoutputPDB, PDBCURinputFile, PDBCURlog)
    os.remove(PDBCURinputFile)
    if not os.path.exists(PDBCURoutputPDB):
        sys.exit('Error 05: Failed to generate Unit Cell PDB file')
    print '\n********** End of Process PDB Section **************************'
    print '****************************************************************'
    print '\n'
    #Parse the newly generated unit cell PDB file
    print '****************************************************************'
    print '********** Parsing PDB Section *********************************\n'
    #return a list of atoms and attributes
    bof, atomList, eof = parsePDB(PDBCURoutputPDB)
    unitCell = getUnitCellParams(pathToPDB)
    print '\n********** End of Parsing PDB Section **************************'
    print '****************************************************************'
    print '\n'
    #Translate the Unit Cell and append the new atom locations to the atomList
    print '****************************************************************'
    print '********** Translate Unit Cell Section *************************\n'
    cartesianVectors = convertToCartesian(unitCell)
    #create shallow copy of the list of atoms to which all translated atomic positions will be added
    transAtomList = copy.copy(atomList) 
    #loop through running the translation subroutine for all combinations of 
    #translations +/- 1 unit cell in a, b and c directions
    for a in range (-1, 2):
        for b in range (-1, 2):
            for c in range (-1, 2):
                #don't translate cells ine the case of an identity translation
                identityTranslation = False
                if a == 0:
                    if b == 0:
                        if c == 0:
                            identityTranslation = True
                if not identityTranslation:
                    #Translate all atoms in the unit cell
                    newTransAtoms = []
                    newTransAtoms = translateUnitCell(atomList, cartesianVectors, a, b, c)
                    #append the translated atom object to list
                    for atm in newTransAtoms:
                        transAtomList.append(copy.copy(atm))
    print ''
    if createAllUnitCellsPDB == 1:
        aucPDBfilepath = '%sAllUnitCells.pdb' % PDBdirectory
        makePDB(bof, transAtomList, eof, aucPDBfilepath)
    print '\n********** End of Translate Unit Cell Section ******************'
    print '****************************************************************'
    print '\n'
    #Discard atoms too far from the asymmetric unit
    print '****************************************************************'
    print '********** Trim Crystal Section ********************************\n'
    #parse a new set of atomic coordinates from the provided asymmetric unit file    
    bof1, auAtomList, eof1 = parsePDB(pathToPDB)
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
    #add PDT to auParams in all dimensions
    print 'Now creating a box surrounding the atoms'
    print 'We only want to consider atoms in this box when we calculate the packing density'
    keepParams = convertParams(auParams, PDT)
    #discard all atoms not in cube defined above
    trimmedAtomList = trimAtoms(transAtomList, keepParams)
    #create PDB file of retained atoms
    if createTrimmedAtomsPDB == 1:
        taPDBfilepath = '%sTrimmedAtoms.pdb' % PDBdirectory
        makePDB(bof, trimmedAtomList, eof, taPDBfilepath)
    print '\n********** End of Trim Crystal Section *************************'
    print '****************************************************************'
    print '\n'
    #Calculate the packing density of each atom
    print '****************************************************************'
    print '********** Calculate Packing Density Section *******************\n'
    minPD, maxPD = calcPDT(auAtomList, trimmedAtomList, PDT)
    noOfGroups = binAtoms(auAtomList, binSize, minPD, maxPD)
    print 'Lowest calculated Packing Density was %s' % minPD
    print 'Highest calculated Packing Density was %s' % maxPD
    print 'Atoms separated into %s bins\n' % noOfGroups
    print '\n********** End of Calculate Packing Density Section ************'
    print '****************************************************************'
    print '\n'
    #Calculate the Bdamage value of each atom
    print '****************************************************************'
    print '********** Calculate Bdamage Section ***************************\n'
    print 'Calculating Bdamage values'
    groupNoAtoms, groupAvBfacs = calcBdam(auAtomList, noOfGroups)
    print 'Writing Bdamage data to output file'
    bDamFileName = '%sBdamage.txt' % fileName
    writeBdam(auAtomList, bDamFileName, groupNoAtoms, groupAvBfacs)
    print '\n********** End of Calculate Bdamage Section ********************'
    print '****************************************************************'
    print '\n'
    #inform the user of the time elapsed while the program was run
    runtime = time.time() - start
    minutes = math.floor(runtime/60)
    seconds = math.fmod(runtime,60)
    if seconds ==1:
        print 'Total time taken for program to run was %02.3f second.\n\n' % seconds
    elif minutes == 0:
        print 'Total time taken for program to run was %02.3f seconds.\n\n' % seconds
    elif minutes == 1:
        print 'Total time taken for program to run was %01.0f minute and %02.3f seconds.\n\n' % (minutes,seconds)   
    else:
        print 'Total time taken for program to run was %01.0f minutes and %02.3f seconds.\n\n' % (minutes,seconds)
#end
CalculateBdamage('2BN3')