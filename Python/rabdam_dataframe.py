

# Script to calculate B-damage for protein atoms
import sys
sys.path.insert(0, './Subroutines')
import time
import math
from CalculateBdamage import rabdam_dataframe

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

# read the filename from the CL input
fileCont = open(sys.argv[1], 'r')
# write the string of the function arguments
functionArgs = fileCont.read()
splitArgs = functionArgs.split(',')
pathToPDBlist = []
for item in splitArgs:
    if '=' not in item:
        pathToPDBlist.append(item.strip())

# remove all whitespace from the inputs
functionArgs = functionArgs.replace(' ', '')
functionArgs = functionArgs.replace('\n', '')
functionArgs = functionArgs.replace('\r', '')
# split the input file into its components
splitArgs = functionArgs.split(',')

# initialise the argument defaults
pdtVal = int(14)
binVal = int(10)
hetatmVal = False
addAtomsList = []
removeAtomsList = []
ucVal = False
aucVal = False
taVal = False

# parse the argument string and assign values to RABDAM inputs
for x in xrange(1, len(splitArgs)):
    if splitArgs[x][0:3] == 'PDT':
        pdtArg = splitArgs[x].split('=')
        pdtVal = float(pdtArg[len(pdtArg)-1])
        if pdtVal == int(pdtVal):
            pdtVal = int(pdtVal)

    elif splitArgs[x][0:7] == 'binSize':
        binArg = splitArgs[x].split('=')
        binVal = int(binArg[len(binArg)-1])

    elif splitArgs[x][0:6] == 'HETATM':
        hetatmArg = splitArgs[x].split('=')
        hetatmVal = hetatmArg[len(hetatmArg) - 1]
        if hetatmVal == 'Keep':
            hetatmVal = True
        if hetatmVal == 'Remove':
            hetatmVal = False

    elif splitArgs[x][0:8] == 'addAtoms':
        addAtomsList = []
        addAtomsArg = splitArgs[x].split('=')
        addAtomsStr = str(addAtomsArg[len(addAtomsArg) - 1])
        if addAtomsStr == '':
            addAtomsList = []
        else:
            addAtomsSubList = addAtomsStr.split(';')
            for item in addAtomsSubList:
                if '-' in item:
                    addAtomsSubList.remove(item)
                    addAtomsRange = item.split('-')
                    if addAtomsRange[-1] != '' and addAtomsRange[0] != '':
                        min_val = int(addAtomsRange[-2])
                        max_val = int(addAtomsRange[-1])
                        addAtomsRange = range(min_val, (max_val + 1))
                        for number in addAtomsRange:
                            addAtomsList.append(str(number))
                    elif addAtomsRange[-1] == '':
                        addAtomsList.append(addAtomsRange[-2])
                    elif addAtomsRange[0] == '':
                        addAtomsList.append(addAtomsRange[-1])
                elif '-' not in item:
                    addAtomsList.append(item)

    elif splitArgs[x][0:11] == 'removeAtoms':
        removeAtomsList = []
        removeAtomsArg = splitArgs[x].split('=')
        removeAtomsStr = str(removeAtomsArg[len(removeAtomsArg) - 1])
        if removeAtomsStr == '':
            removeAtomsList = []
        else:
            removeAtomsSubList = removeAtomsStr.split(';')
            for item in removeAtomsSubList:
                if '-' in item:
                    removeAtomsSubList.remove(item)
                    removeAtomsRange = item.split('-')
                    if removeAtomsRange[-1] != '' and removeAtomsRange[0] != '':
                        min_val = int(removeAtomsRange[-2])
                        max_val = int(removeAtomsRange[-1])
                        removeAtomsRange = range(min_val, (max_val + 1))
                        for number in removeAtomsRange:
                            removeAtomsList.append(str(number))
                    elif removeAtomsRange[-1] == '':
                        removeAtomsList.append(removeAtomsRange[-2])
                    elif removeAtomsRange[0] == '':
                        removeAtomsList.append(removeAtomsRange[-1])
                elif '-' not in item:
                    removeAtomsList.append(item)

    elif splitArgs[x][0:11] == 'createAUpdb':
        auArg = splitArgs[x].split('=')
        auVal = auArg[len(auArg) - 1]
        if auVal == 'True':
            auVal = True
        elif auVal == 'False':
            auVal = False
        else:
            sys.exit('Error 00: Input file is formatted incorrectly\nRead the handbook and amend the input file\n(Try looking at the createAUpdb argument)')

    elif splitArgs[x][0:11] == 'createUCpdb':
        ucArg = splitArgs[x].split('=')
        ucVal = ucArg[len(ucArg) - 1]
        if ucVal == 'True':
            ucVal = True
        elif ucVal == 'False':
            ucVal = False
        else:
            sys.exit('Error 00: Input file is formatted incorrectly\nRead the handbook and amend the input file\n(Try looking at the createUCpdb argument)')

    elif splitArgs[x][0:12] == 'createAUCpdb':
        aucArg = splitArgs[x].split('=')
        aucVal = aucArg[len(aucArg) - 1]
        if aucVal == 'True':
            aucVal = True
        elif aucVal == 'False':
            aucVal = False
        else:
            sys.exit('Error 00: Input file is formatted incorrectly\nRead the handbook and amend the input file\n(Try looking at the createAUCpdb argument)')

    elif splitArgs[x][0:11] == 'createTApdb':
        taArg = splitArgs[x].split("=")
        taVal = taArg[len(taArg) - 1]
        if taVal == 'True':
            taVal = True
        elif taVal == 'False':
            taVal = False
        else:
            sys.exit('Error 00: Input file is formatted incorrectly\nRead the handbook and amend the input file\n(Try looking at the createTApdb argument)')

# run function with arguments from input file
for item in pathToPDBlist:
    rabdam_dataframe(item, PDT=pdtVal, binSize=binVal, HETATM=hetatmVal, addAtoms=addAtomsList, removeAtoms=removeAtomsList, createAUpdb=auVal, createUCpdb=ucVal, createAUCpdb=aucVal, createTApdb=taVal)

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
