# Script to calculate B-damage for protein atom
# Copyright 2015 Thomas Dixon
print('Copyright 2015 Thomas Dixon\n')
print('\n')
#import packages required for running the program
import time
#Input: the file path to the pdb for which you want to calculate B-damage factors, the 'Packing Density Threshold' (Angstroms) and bin size

def CalculateBdamage(pathToPDB, PDT=14, binSize=10):
    start = time.time()
    print '################################################################\n'
    print '################################################################\n'
    print '################################################################\n'
    print '### Program to calculate B Damage                            ###\n'
    print '################################################################\n'
    print '\n\n'
    startIndex = time.gmtime()
    year = startIndex.tm_year
    month = startIndex.tm_mon
    day = startIndex.tm_mday
    hour = startIndex.tm_hour
    minute = startIndex.tm_min
    second = startIndex.tm_sec
    #inform the user of the start time of the program
    print 'This program was run on %f\/%f\/%f at %f\:%f\:%f\n' % day % month % year % hour % minute % second
    
    #inform the user of the time elapsed while the program was run
    runtime = time.time() - start
    print 'Program run in %f seconds\n' % runtime
#    print 'Total time taken for program to run was %.0f minutes and %.0f seconds.\n\n' % floor(timetaken/60) % rem(timetaken/60)
    
CalculateBdamage('2CCC')
