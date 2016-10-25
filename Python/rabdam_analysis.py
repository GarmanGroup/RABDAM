

import sys
import time
import math
import copy

sys.path.insert(0, './Subroutines')
duplicate = copy.copy

from CalculateBdamage import rabdam_analysis

# An outer layer to the pipeline scripts. This allows the B_damage analysis
# pipeline to be run from the command line by calling:
#
#          python rabdam_analysis.py INPUT.txt

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
pathToPDBlist = []
for item in splitArgs:
    if '=' not in item:
        pathToPDBlist.append(item.strip())
pathToPDBlist = filter(None, pathToPDBlist)

# Reads in the remaining rabdam_analysis function arguments from INPUT.txt.
functionArgs = functionArgs.replace(' ', '')
functionArgs = functionArgs.replace('\n', '')
functionArgs = functionArgs.replace('\r', '')
splitArgs = functionArgs.split(',')
splitArgsDuplicate = duplicate(splitArgs)
for item in splitArgsDuplicate:
    if '=' not in item:
        splitArgs.remove(item)

# Initialises argument default values.
thresholdVal = float(0.02)
highlightAtomsList = []
run = 'rabdam_analysis'

# Assigns argument values as provided in INPUT.txt.
for x in xrange(0, len(splitArgs)):
    if splitArgs[x][0:9] == 'threshold':
        thresholdArg = splitArgs[x].split('=')
        thresholdVal = thresholdArg[len(thresholdArg)-1]
        if '%' in thresholdVal:
            thresholdVal = thresholdVal.replace('%', '')
            thresholdVal = float(thresholdVal) / 100
        thresholdVal = float(thresholdVal)

    elif splitArgs[x][0:14] == 'highlightAtoms':
        highlightAtomsList = []
        highlightAtomsArg = splitArgs[x].split('=')
        highlightAtomsStr = str(highlightAtomsArg[len(highlightAtomsArg)-1])
        if highlightAtomsStr == '':
            highlightAtomsList = []
        else:
            highlightAtomsSubList = highlightAtomsStr.split(';')
            for item in highlightAtomsSubList:
                highlightAtomsRange = str(item).split('-')
                if len(highlightAtomsRange) == 2:
                    min_val = int(highlightAtomsRange[-2])
                    max_val = int(highlightAtomsRange[-1])
                    highlightAtomsRange = range(min_val, (max_val+1))
                    for number in highlightAtomsRange:
                        highlightAtomsList.append(str(number))
                elif len(highlightAtomsRange) == 1:
                    highlightAtomsList.append(highlightAtomsRange[-1])

for item in pathToPDBlist:
    # Generates output analysis files from pre-calculated B_damage values.
    rabdam_analysis(
        item, threshold=thresholdVal, highlightAtoms=highlightAtomsList,
        run=run
        )

runtime = time.time() - start
mins = math.floor(runtime/60)
secs = math.fmod(runtime, 60)
if mins == 0:
    if secs == 1:
        print 'Program run time = %02.3f second\n\n' % secs
    else:
        print 'Programe run time = %02.3f seconds\n\n' % secs
elif mins == 1:
    if secs == 1:
        print 'Program run time = %01.0f minute, %02.3f second\n\n' % (mins,
                                                                       secs)
    else:
        print 'Program run time = %01.0f minute, %02.3f seconds\n\n' % (mins,
                                                                        secs)
else:
    print 'Program run time = %01.0f minutes, %02.3f seconds\n\n' % (mins,
                                                                     secs)
