

# Script to calculate B-damage for protein atoms
import sys
sys.path.insert(0, './Subroutines')
from CalculateBdamage import rabdam

# read the filename from the CL input
fileCont = open(sys.argv[1], 'r')
# write the string of the function arguments
functionArgs = fileCont.read()
# remove all whitespace from the inputs
functionArgs = functionArgs.replace(' ', '')
functionArgs = functionArgs.replace('\n', '')
functionArgs = functionArgs.replace('\r', '')
# split the input file into its components
splitArgs = functionArgs.split(',')
# parse the argument string and assign values to RABDAM inputs
pathToPDB = splitArgs[0]
# initialise the argument defaults
pdtVal = int(14)
binVal = int(10)
addAtomsList = []
removeAtomsList = []
thresholdVal = float(0.02)
aucVal = False
taVal = False

for x in xrange(1, len(splitArgs)):
    if splitArgs[x][0:3] == 'PDT':
        pdtArg = splitArgs[x].split('=')
        pdtVal = float(pdtArg[len(pdtArg)-1])
        if pdtVal == int(pdtVal):
            pdtVal = int(pdtVal)
        break
    elif splitArgs[x][0:7] == 'binSize':
        binArg = splitArgs[x].split('=')
        binVal = int(binArg[len(binArg)-1])
        break
    elif splitArgs[x][0:8] == 'addAtoms':
        addAtomsArg = splitArgs[x].split('=')
        addAtomsStr = str(addAtomsArg[len(addAtomsArg) - 1])
        addAtomsList = addAtomsStr.split(';')
        break
    elif splitArgs[x][0:11] == 'removeAtoms':
        removeAtomsArg = splitArgs[x].split('=')
        removeAtomsStr = str(removeAtomsArg[len(removeAtomsArg) - 1])
        removeAtomsList = removeAtomsStr.split(';')
        break
    elif splitArgs[x][0:9] == 'threshold':
        thresholdArg = splitArgs[x].split('=')
        thresholdVal = float(thresholdArg[len(thresholdArg) - 1])
        break
    elif splitArgs[x][0:12] == 'createAUCpdb':
        aucArg = splitArgs[x].split('=')
        aucVal = aucArg[len(aucArg)-1]
        if aucVal == 'True':
            aucVal = True
        elif aucVal == 'False':
            aucVal = False
        else:
            sys.exit('Error 00: Input file is formatted incorrectly\nRead the handbook and amend the input file\n(Try looking at the createAUCpdb argument)')
        break
    elif splitArgs[x][0:11] == 'createTApdb':
        taArg = splitArgs[x].split("=")
        taVal = taArg[(len(taArg))-1]
        if taVal == 'True':
            taVal = True
        elif taVal == 'False':
            taVal = False
        else:
            sys.exit('Error 00: Input file is formatted incorrectly\nRead the handbook and amend the input file\n(Try looking at the createTApdb argument)')
        break
    else:
        sys.exit('Error 00: Input file is formatted incorrectly\nRead the handbook and amend the input file')

# run fucntion with arguments from input file
rabdam(pathToPDB, PDT=pdtVal, binSize=binVal, addAtoms=addAtomsList, removeAtoms=removeAtomsList, threshold=thresholdVal, createAUCpdb=aucVal, createTApdb=taVal)
