

import sys
import time
import math

sys.path.insert(0, './Subroutines')

from CalculateBdamage import rabdam_dataframe

# An outer layer to the pipeline scripts. This allows the B_damage calculation
# pipeline to be run from the command line by calling:
#
#          python rabdam_dataframe.py INPUT.txt

start = time.time()
startIndex = time.gmtime()
year = startIndex.tm_year
month = startIndex.tm_mon
day = startIndex.tm_mday
hour = startIndex.tm_hour
minute = startIndex.tm_min
second = startIndex.tm_sec

print 'This program was run on %d/%d/%d at %02.0f:%02.0f:%02.0f\n\n' % (
    day, month, year, hour, minute, second
    )

# Reads in the PDB file names listed in INPUT.txt.
fileCont = open(sys.argv[1], 'r')
functionArgs = fileCont.read()
splitArgs = functionArgs.split(',')
pathToPDBlist = [item.strip() for item in splitArgs if '=' not in item]
pathToPDBlist = filter(None, pathToPDBlist)
if len(pathToPDBlist) == 0:
    sys.exit('No input PDB code / file listed in INPUT.txt')

# Reads in the remaining rabdam_dataframe function arguments from INPUT.txt.
functionArgs = functionArgs.replace(' ', '')
functionArgs = functionArgs.replace('\n', '')
functionArgs = functionArgs.replace('\r', '')
splitArgs = functionArgs.split(',')
splitArgs = [item for item in splitArgs if '=' in item]

# Initialises argument default values.
pdtVal = int(14)
windowVal = float(0.02)
protOrNAVal = 'BOTH'
hetatmVal = False
addAtomsList = []
removeAtomsList = []
auVal = False
ucVal = False
aucVal = False
taVal = False
run = 'rabdam_dataframe'

# Assigns argument values as provided in INPUT.txt.
for x in xrange(0, len(splitArgs)):
    if splitArgs[x][0:3] == 'PDT':
        pdtArg = splitArgs[x].split('=')
        pdtVal = float(pdtArg[len(pdtArg)-1])

    elif splitArgs[x][0:10] == 'windowSize':
        windowArg = splitArgs[x].split('=')
        windowVal = windowArg[len(windowArg)-1]
        if '%' in windowVal:
            windowVal = windowVal.replace('%', '')
            windowVal = float(windowVal) / 100
        windowVal = float(windowVal)

    elif splitArgs[x][0:20] == 'proteinOrNucleicAcid':
        protOrNAArg = splitArgs[x].split('=')
        protOrNAVal = str((protOrNAArg[len(protOrNAArg)-1]).upper())

    elif splitArgs[x][0:6] == 'HETATM':
        hetatmArg = splitArgs[x].split('=')
        hetatmVal = str(hetatmArg[len(hetatmArg)-1].upper())
        if hetatmVal == 'KEEP':
            hetatmVal = True
        elif hetatmVal == 'REMOVE':
            hetatmVal = False

    elif splitArgs[x][0:8] == 'addAtoms':
        addAtomsList = []
        addAtomsArg = splitArgs[x].split('=')
        addAtomsStr = str((addAtomsArg[len(addAtomsArg)-1]).upper())
        if addAtomsStr != '':
            addAtomsSubList = addAtomsStr.split(';')
            for item in addAtomsSubList:
                addAtomsRange = str(item).split('-')
                if len(addAtomsRange) == 2:
                    min_val = int(addAtomsRange[-2])
                    max_val = int(addAtomsRange[-1])
                    addAtomsRange = range(min_val, (max_val+1))
                    for number in addAtomsRange:
                        addAtomsList.append(str(number))
                elif len(addAtomsRange) == 1:
                    addAtomsList.append(addAtomsRange[-1])

    elif splitArgs[x][0:11] == 'removeAtoms':
        removeAtomsList = []
        removeAtomsArg = splitArgs[x].split('=')
        removeAtomsStr = str((removeAtomsArg[len(removeAtomsArg)-1]).upper())
        if removeAtomsStr != '':
            removeAtomsSubList = removeAtomsStr.split(';')
            for item in removeAtomsSubList:
                removeAtomsRange = str(item).split('-')
                if len(removeAtomsRange) == 2:
                    min_val = int(removeAtomsRange[-2])
                    max_val = int(removeAtomsRange[-1])
                    removeAtomsRange = range(min_val, (max_val+1))
                    for number in removeAtomsRange:
                        removeAtomsList.append(str(number))
                elif len(removeAtomsRange) == 1:
                    removeAtomsList.append(removeAtomsRange[-1])

    elif splitArgs[x][0:11] == 'createAUpdb':
        auArg = splitArgs[x].split('=')
        auVal = auArg[len(auArg)-1].upper()
        if auVal == 'TRUE':
            auVal = True
        elif auVal == 'FALSE':
            auVal = False

    elif splitArgs[x][0:11] == 'createUCpdb':
        ucArg = splitArgs[x].split('=')
        ucVal = ucArg[len(ucArg)-1].upper()
        if ucVal == 'TRUE':
            ucVal = True
        elif ucVal == 'FALSE':
            ucVal = False

    elif splitArgs[x][0:12] == 'createAUCpdb':
        aucArg = splitArgs[x].split('=')
        aucVal = aucArg[len(aucArg)-1].upper()
        if aucVal == 'TRUE':
            aucVal = True
        elif aucVal == 'FALSE':
            aucVal = False

    elif splitArgs[x][0:11] == 'createTApdb':
        taArg = splitArgs[x].split('=')
        taVal = taArg[len(taArg)-1].upper()
        if taVal == 'TRUE':
            taVal = True
        elif taVal == 'FALSE':
            taVal = False

for item in pathToPDBlist:
    # Calculates B_damage values and writes them to a DataFrame.
    rabdam_dataframe(
        item, PDT=pdtVal, windowSize=windowVal, protOrNA=protOrNAVal,
        HETATM=hetatmVal, addAtoms=addAtomsList, removeAtoms=removeAtomsList,
        createAUpdb=auVal, createUCpdb=ucVal, createAUCpdb=aucVal,
        createTApdb=taVal, run=run
        )

runtime = time.time() - start
mins = math.floor(runtime/60)
secs = math.fmod(runtime, 60)
if mins == 0:
    print 'Program run time: %02.3f sec\n\n' % secs
elif mins >= 1:
    print 'Program run time: %01.0f min, %02.3f sec\n\n' % (mins, secs)
