def CalculateBdamage(pathToPDB, PDT=14, binSize=10, createUnitCellPDB=0, createTrimmedUnitCellPDB=1):
    # Script to calculate B-damage for protein atom
    # Copyright 2015 Thomas Dixon
    print('\n')
    print('Copyright 2015 Thomas Dixon\n')
    #import packages required for running the program
    import time
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
    print 'This program was run on %d/%d/%d at %02.0f:%02.0f:%02.0f\n\n'%(day, month, year, hour, minute, second)
    print '################################################################'
    print '################################################################'
    print '### Program to calculate B Damage ##############################'
    print '################################################################'
    print '\n'
    #inform the user of the inputs used in the program
    print '****************************************************************'
    print '************************ Input Section *************************\n'
    print 'Calculating B Damage for %s' % pathToPDB
    if PDT == 14:
        print 'No packing density threshold supplied, using default value of 14 Angstroms'
    else:
        print 'Packing density threshold defined by user at %s Angstroms' % PDT
    if binSize == 10:
        print 'No bin size supplied, using default value of 10'
    else:
        print 'Bin size defined by user at %s\n' % binSize
    print '******************* End of Input Section ***********************'
    print '****************************************************************'
    print '\n'
    #Process pdb file with PDBCUR to generate a unit cell and then store the name of the processed pdb
    print '****************************************************************'
    print '********************* Process PDB Section **********************\n'
    print 'Processing PDB file to:'
    print '1 - Remove hydrogen atoms.'
    print '2 - Keep most probable alternate conformation.'
    print '3 - Remove anisotropic records from the file.'
    print '4 - Generate the Unit Cell from Symmetry Operations.\n'
    time.sleep(2)
    #inform the user of the time elapsed while the program was run
    runtime = time.time() - start
    print 'Program run in %f seconds\n' % runtime
#    print 'Total time taken for program to run was %.0f minutes and %.0f seconds.\n\n' % floor(timetaken/60) % rem(timetaken/60)    

#end