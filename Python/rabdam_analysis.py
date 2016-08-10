

# Script to calculate B-damage for protein atoms
import sys
sys.path.insert(0, './Subroutines')
import time
import math
from CalculateBdamage import rabdam_analysis

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
thresholdVal = float(0.02)
highlightAtomsList = []

# parse the argument string and assign values to RABDAM inputs
for x in xrange(1, len(splitArgs)):
    if splitArgs[x][0:9] == 'threshold':
        thresholdArg = splitArgs[x].split('=')
        thresholdVal = float(thresholdArg[len(thresholdArg) - 1])

    elif splitArgs[x][0:14] == 'highlightAtoms':
        highlightAtomsArg = splitArgs[x].split('=')
        highlightAtomsStr = str(highlightAtomsArg[len(highlightAtomsArg) - 1])
        highlightAtomsList = highlightAtomsStr.split(';')

# run fucntion with arguments from input file
for item in pathToPDBlist:
    rabdam_analysis(item, threshold=thresholdVal, highlightAtoms=highlightAtomsList)

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
