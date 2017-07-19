

import sys
import os
import argparse
import time
import math
import numpy as np

sys.path.insert(0, './Subroutines')

from CalculateBdamage import rabdam

# An outer layer to the pipeline scripts. This allows the complete pipeline
# to be run from the command line by calling:
#
#          python rabdam.py INPUT.txt

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

# Reads in command line inputs
parser = argparse.ArgumentParser()
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument('-i', '--input', help='Path to input file specifying '
                   'program options')
group.add_argument('-f', '--file_name', nargs='+', help='Specifies input pdb '
                   'file for Bdamage analysis - this option allows the RABDAM '
                   'program to be run (using default parameter values) '
                   'without providing an input file specifying program '
                   'options')
parser.add_argument('-o', '--output', help='Specifies whether to run the '
                    'complete program (default), to calculate Bdamage values '
                    'only ("df" / "dataframe"), or to analyse pre-calculated '
                    'Bdamage values only ("analysis")')
args = parser.parse_args()

if vars(args)['input'] is not None:
    cwd = os.getcwd()
    # Reads in the PDB file name(s) listed in INPUT.txt.
    input_file_loc = vars(args)['input']
    input_file_loc = input_file_loc.replace('\'', '')
    input_file_loc_list = input_file_loc.split('/')
    input_file = input_file_loc_list[len(input_file_loc_list)-1]
    if len(input_file_loc_list) > 1:
        input_file_loc = input_file_loc.replace(input_file, '')
    else:
        input_file_loc = cwd
    os.chdir(input_file_loc)
    fileCont = open(input_file, 'r')
    functionArgs = fileCont.read()
    fileCont.close()
    os.chdir(cwd)
    splitArgs = functionArgs.split(',')
    pathToPDBlist = [item.strip() for item in splitArgs if '=' not in item]
    pathToPDBlist = filter(None, pathToPDBlist)

    # Reads in the remaining rabdam function arguments from INPUT.txt.
    functionArgs = functionArgs.replace(' ', '')
    functionArgs = functionArgs.replace('\n', '')
    functionArgs = functionArgs.replace('\r', '')
    splitArgs = functionArgs.split(',')
    splitArgs = [item for item in splitArgs if '=' in item]

elif vars(args)['file_name'] is not None:
    # Reads in PDB file name(s) listed on command line input
    pathToPDBlist = vars(args)['file_name']
    splitArgs = []

if len(pathToPDBlist) == 0:
    sys.exit('No input PDB code / file provided')

# Initialises argument default values.
outputLoc = '.'
pdtList = [14]
windowList = [0.02]
protOrNAVal = 'BOTH'
hetatmVal = False
addAtomsList = []
removeAtomsList = []
thresholdVal = float(0.02)
highlightAtomsList = []
auVal = False
ucVal = False
aucVal = False
taVal = False
run = 'rabdam'

# Assigns argument values as provided in INPUT.txt.
for x in xrange(0, len(splitArgs)):
    if splitArgs[x][0:3].lower() == 'dir':
        dirArg = splitArgs[x].split('=')
        outputLoc = dirArg[len(dirArg)-1]
        if outputLoc == '':
            print 'WAAAAAAAA!'
            outputLoc = '.'
        outputLoc = outputLoc.replace('\'', '')

    elif splitArgs[x][0:3].lower() == 'pdt':
        pdtArg = splitArgs[x].split('=')
        pdtArg = pdtArg[len(pdtArg)-1]
        if '-' in pdtArg:
            pdtRange = pdtArg.split('-')
            min_pdt = float(pdtRange[0])
            max_pdt = float(pdtRange[len(pdtRange)-1])
            pdtList = np.arange(min_pdt, max_pdt, 0.5)
        else:
            pdtList.append(float(pdtArg))

    elif splitArgs[x][0:10].lower() == 'windowsize':
        windowArg = splitArgs[x].split('=')
        windowArg = windowArg[len(windowArg)-1]
        windowVals = windowArg.split(';')
        windowList = []
        for number in windowVals:
            if '%' in number:
                number = number.replace('%', '')
                number = float(number) / 100
                windowList.append(number)
            else:
                windowList.append(float(number))

    elif splitArgs[x][0:20].lower() == 'proteinornucleicacid':
        protOrNAArg = splitArgs[x].split('=')
        protOrNAVal = str((protOrNAArg[len(protOrNAArg)-1]).upper())

    elif splitArgs[x][0:6].lower() == 'hetatm':
        hetatmArg = splitArgs[x].split('=')
        hetatmVal = str((hetatmArg[len(hetatmArg)-1]).lower())
        if hetatmVal == 'keep':
            hetatmVal = True
        elif hetatmVal == 'remove':
            hetatmVal = False

    elif splitArgs[x][0:8].lower() == 'addatoms':
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

    elif splitArgs[x][0:11].lower() == 'removeatoms':
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

    elif splitArgs[x][0:9].lower() == 'threshold':
        thresholdArg = splitArgs[x].split('=')
        thresholdVal = thresholdArg[len(thresholdArg)-1]
        if '%' in thresholdVal:
            thresholdVal = thresholdVal.replace('%', '')
            thresholdVal = float(thresholdVal) / 100
        thresholdVal = float(thresholdVal)

    elif splitArgs[x][0:14].lower() == 'highlightatoms':
        highlightAtomsList = []
        highlightAtomsArg = splitArgs[x].split('=')
        highlightAtomsStr = str(highlightAtomsArg[len(highlightAtomsArg)-1])
        if highlightAtomsStr != '':
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

    elif splitArgs[x][0:11].lower() == 'createaupdb':
        auArg = splitArgs[x].split('=')
        auVal = auArg[len(auArg)-1].lower()
        if auVal == 'true':
            auVal = True
        elif auVal == 'false':
            auVal = False

    elif splitArgs[x][0:11].lower() == 'createucpdb':
        ucArg = splitArgs[x].split('=')
        ucVal = ucArg[len(ucArg)-1].lower()
        if ucVal == 'true':
            ucVal = True
        elif ucVal == 'false':
            ucVal = False

    elif splitArgs[x][0:12].lower() == 'createaucpdb':
        aucArg = splitArgs[x].split('=')
        aucVal = aucArg[len(aucArg)-1].lower()
        if aucVal == 'true':
            aucVal = True
        elif aucVal == 'false':
            aucVal = False

    elif splitArgs[x][0:11].lower() == 'createtapdb':
        taArg = splitArgs[x].split('=')
        taVal = taArg[len(taArg)-1].lower()
        if taVal == 'true':
            taVal = True
        elif taVal == 'false':
            taVal = False

for item in pathToPDBlist:
    for windowVal in windowList:
        for pdtVal in pdtList:
            y = rabdam(pathToPDB=item, outputDir=outputLoc, PDT=pdtVal,
                       windowSize=windowVal, protOrNA=protOrNAVal,
                       HETATM=hetatmVal, addAtoms=addAtomsList,
                       removeAtoms=removeAtomsList, threshold=thresholdVal,
                       highlightAtoms=highlightAtomsList,
                       createAUpdb=auVal, createUCpdb=ucVal,
                       createAUCpdb=aucVal, createTApdb=taVal)
            if vars(args)['output'] is None:
                # Calculates B_damage values and writes them to a DataFrame.
                y.rabdam_dataframe(run='rabdam')
                # Generates output analysis files from pre-calculated B_damage
                # values.
                y.rabdam_analysis(run='rabdam')
            elif vars(args)['output'].lower() in ['dataframe', 'df']:
                # Calculates B_damage values and writes them to a DataFrame.
                y.rabdam_dataframe(run='rabdam_dataframe')
            elif vars(args)['output'].lower() in ['analysis']:
                # Generates output analysis files from pre-calculated B_damage
                # values.
                y.rabdam_analysis(run='rabdam_analysis')


runtime = time.time() - start
mins = math.floor(runtime/60)
secs = math.fmod(runtime, 60)
if mins == 0:
    print 'Program run time: %02.3f sec\n\n' % secs
elif mins >= 1:
    print 'Program run time: %01.0f min %02.3f sec\n\n' % (mins, secs)
