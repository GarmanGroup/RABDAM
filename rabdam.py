

# An outer layer to the pipeline scripts. Depending upon the flags specified
# in the command line input, this script will run either the complete / a
# subsection of the pipeline.

import sys
import os
import argparse
import time
import math
import numpy as np

sys.path.insert(0, './Subroutines')

from CalculateBdamage import rabdam


# Initialises program start time
start = time.time()
startIndex = time.localtime()
year = startIndex.tm_year
month = startIndex.tm_mon
day = startIndex.tm_mday
hour = startIndex.tm_hour
minute = startIndex.tm_min
second = startIndex.tm_sec

print '\nThis program was run on %d/%d/%d at %02.0f:%02.0f:%02.0f\n\n' % (
    day, month, year, hour, minute, second
    )

# Reads in command line inputs. There are three recognised flags: -i, -f and
# -o. Program inputs are specified by either the -i or the -f flag; provision
# of one of these flags is compulsory. Program outputs are specified by the
# optional -o flag.
parser = argparse.ArgumentParser()
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument('-i', '--input', help='Path to input file listing program '
                   'options')
group.add_argument('-f', '--pdb_file', nargs='+', help='Specifies input pdb '
                   'file for BDamage analysis - this option allows the RABDAM '
                   'program to be run (using default parameter values) '
                   'without providing an input file listing program options')
parser.add_argument('-o', '--output', help='Specifies whether to run the '
                    'complete program (default), to calculate BDamage values '
                    'only ("df" / "dataframe"), or to analyse pre-calculated '
                    'BDamage values only ("analysis")')
args = parser.parse_args()

cwd = os.getcwd()
# Reads in program options from input file specified by -i flag
if vars(args)['input'] is not None:
    # Reads in PDB file name(s) listed in input file
    input_file_loc = vars(args)['input']
    input_file_loc = input_file_loc.replace('\\', '/')
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
    pathToPDBlist = [item.strip() for item in splitArgs if '=' not in item
                     and item is not None]
    # Reads in remaining program options specified in input file
    functionArgs = functionArgs.replace(' ', '')
    functionArgs = functionArgs.replace('\n', '')
    functionArgs = functionArgs.replace('\r', '')
    splitArgs = functionArgs.split(',')
    splitArgs = [item for item in splitArgs if '=' in item]

# Reads in PDB file name(s) specified by -f flag listed in command line input
elif vars(args)['pdb_file'] is not None:
    pathToPDBlist = vars(args)['pdb_file']
    splitArgs = []

if len(pathToPDBlist) == 0:
    sys.exit('No input PDB code / file provided')

# Initialises default program options
outputLoc = cwd
pdtList = [14.0]
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

# If input file is provided, program options are updated to the values it
# specifies
for x in xrange(0, len(splitArgs)):
    # Specifies location to which output 'Logfiles' directory is written
    if splitArgs[x][0:3].lower() == 'dir':
        dirArg = splitArgs[x].split('=')
        outputLoc = dirArg[len(dirArg)-1]
        if outputLoc == '':
            outputLoc = cwd
        outputLoc = outputLoc.replace('\\', '/')

    # Specifies packing density threshold
    elif splitArgs[x][0:3].lower() == 'pdt':
        pdtArg = splitArgs[x].split('=')
        pdtArg = pdtArg[len(pdtArg)-1]
        if '-' in pdtArg:
            pdtRange = pdtArg.split('-')
            min_pdt = float(pdtRange[0])
            max_pdt = float(pdtRange[len(pdtRange)-1])
            pdtList = np.arange(min_pdt, max_pdt, 0.5)
        else:
            pdtList = [float(pdtArg)]

    # Specifies size of sliding window (as a percentage of the total number of
    # atoms considered for BDamage analysis)
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

    # Specifies whether to include protein atoms alone, nucleic atoms alone,
    # or both atom types (if present) in the BDamage calculation
    elif splitArgs[x][0:20].lower() == 'proteinornucleicacid':
        protOrNAArg = splitArgs[x].split('=')
        protOrNAVal = protOrNAArg[len(protOrNAArg)-1].upper()

    # Specifies whether to remove HETATM from the BDamage calculation or to
    # retain them
    elif splitArgs[x][0:6].lower() == 'hetatm':
        hetatmArg = splitArgs[x].split('=')
        hetatmVal = hetatmArg[len(hetatmArg)-1].lower()
        if hetatmVal == 'keep':
            hetatmVal = True
        elif hetatmVal == 'remove':
            hetatmVal = False

    # Lists atoms (via either their atom numbers or their residue names) to be
    # included in the BDamage calculation. (This is useful to add back in a
    # subset of atoms of interest that has been removed by the
    # proteinOrNucleicAcid / HETATM / removeAtoms options.)
    elif splitArgs[x][0:8].lower() == 'addatoms':
        addAtomsList = []
        addAtomsArg = splitArgs[x].split('=')
        addAtomsStr = addAtomsArg[len(addAtomsArg)-1].upper()
        if addAtomsStr != '':
            addAtomsSubList = addAtomsStr.split(';')
            for item in addAtomsSubList:
                addAtomsRange = item.split('-')
                if len(addAtomsRange) == 2:
                    min_val = int(addAtomsRange[-2])
                    max_val = int(addAtomsRange[-1])
                    addAtomsRange = xrange(min_val, (max_val+1))
                    for number in addAtomsRange:
                        addAtomsList.append(str(number))
                elif len(addAtomsRange) == 1:
                    addAtomsList.append(addAtomsRange[-1])

    # Lists atoms (via either their atom numbers or their residue names) to be
    # removed from the BDamage calculation. (This is useful to remove
    # additional atoms not covered by the proteinOrNucleicAcid / HETATM
    # options.)
    elif splitArgs[x][0:11].lower() == 'removeatoms':
        removeAtomsList = []
        removeAtomsArg = splitArgs[x].split('=')
        removeAtomsStr = removeAtomsArg[len(removeAtomsArg)-1].upper()
        if removeAtomsStr != '':
            removeAtomsSubList = removeAtomsStr.split(';')
            for item in removeAtomsSubList:
                removeAtomsRange = item.split('-')
                if len(removeAtomsRange) == 2:
                    min_val = int(removeAtomsRange[-2])
                    max_val = int(removeAtomsRange[-1])
                    removeAtomsRange = xrange(min_val, (max_val+1))
                    for number in removeAtomsRange:
                        removeAtomsList.append(str(number))
                elif len(removeAtomsRange) == 1:
                    removeAtomsList.append(removeAtomsRange[-1])

    # Specifies the number of atoms (as a percentage of the total number of
    # atoms considered for BDamage analysis) with the highest BDamage values
    # to be listed in the program output
    elif splitArgs[x][0:9].lower() == 'threshold':
        thresholdArg = splitArgs[x].split('=')
        thresholdVal = thresholdArg[len(thresholdArg)-1]
        if '%' in thresholdVal:
            thresholdVal = thresholdVal.replace('%', '')
            thresholdVal = float(thresholdVal) / 100
        thresholdVal = float(thresholdVal)

    # Lists atoms (via their atom numbers) whose BDamage values are to be
    # indicated on the kernel density estimate of the BDamage distribution
    # output by the program. Note that it is recommended no more than 6 atoms
    # are listed in the highlightAtoms option in the input file (beyond 6
    # atoms, the colour scheme will repeat itself, and in addition the key may
    # not fit onto the graph).
    elif splitArgs[x][0:14].lower() == 'highlightatoms':
        highlightAtomsList = []
        highlightAtomsArg = splitArgs[x].split('=')
        highlightAtomsStr = highlightAtomsArg[len(highlightAtomsArg)-1]
        if highlightAtomsStr != '':
            highlightAtomsSubList = highlightAtomsStr.split(';')
            for item in highlightAtomsSubList:
                highlightAtomsRange = item.split('-')
                if len(highlightAtomsRange) == 2:
                    min_val = int(highlightAtomsRange[-2])
                    max_val = int(highlightAtomsRange[-1])
                    highlightAtomsRange = xrange(min_val, (max_val+1))
                    for number in highlightAtomsRange:
                        highlightAtomsList.append(str(number))
                elif len(highlightAtomsRange) == 1:
                    highlightAtomsList.append(highlightAtomsRange[-1])

    # Specifies whether to create a PDB file of the (filtered) asymmetric unit
    elif splitArgs[x][0:11].lower() == 'createaupdb':
        auArg = splitArgs[x].split('=')
        auVal = auArg[len(auArg)-1].lower()
        if auVal in ['true', 'yes', 't', 'y']:
            auVal = True
        elif auVal in ['false', 'no', 'f', 'n']:
            auVal = False

    # Specifies whether to create a PDB file of the unit cell
    elif splitArgs[x][0:11].lower() == 'createucpdb':
        ucArg = splitArgs[x].split('=')
        ucVal = ucArg[len(ucArg)-1].lower()
        if ucVal in ['true', 'yes', 't', 'y']:
            ucVal = True
        elif ucVal in ['false', 'no', 'f', 'n']:
            ucVal = False

    # Specifies whether to create a PDB file of the 3x3x3 unit cell assembly
    elif splitArgs[x][0:12].lower() == 'createaucpdb':
        aucArg = splitArgs[x].split('=')
        aucVal = aucArg[len(aucArg)-1].lower()
        if aucVal in ['true', 'yes', 't', 'y']:
            aucVal = True
        elif aucVal in ['false', 'no', 'f', 'n']:
            aucVal = False

    # Specifies whether to create a PDB file of the trimmed atoms assembly
    # (which consists of the asymmetric unit, plus every atom in the 3x3x3 unit
    # cell assembly within a PDT (default=14) Angstrom radius of the asymmetric
    # unit)
    elif splitArgs[x][0:11].lower() == 'createtapdb':
        taArg = splitArgs[x].split('=')
        taVal = taArg[len(taArg)-1].lower()
        if taVal in ['true', 'yes', 't', 'y']:
            taVal = True
        elif taVal in ['false', 'no', 'f', 'n']:
            taVal = False

# Runs the BDamage calculation for every specified PDB file
for item in pathToPDBlist:
    for windowVal in windowList:
        for pdtVal in pdtList:
            # Initialises rabdam object
            pdb = rabdam(pathToPDB=item, outputDir=outputLoc, PDT=pdtVal,
                         windowSize=windowVal, protOrNA=protOrNAVal,
                         HETATM=hetatmVal, addAtoms=addAtomsList,
                         removeAtoms=removeAtomsList, threshold=thresholdVal,
                         highlightAtoms=highlightAtomsList,
                         createAUpdb=auVal, createUCpdb=ucVal,
                         createAUCpdb=aucVal, createTApdb=taVal)
            if vars(args)['output'] is None:
                # Runs full program
                pdb.rabdam_dataframe(run='rabdam')
                pdb.rabdam_analysis(run='rabdam')
            elif vars(args)['output'].lower() in ['dataframe', 'df']:
                # Runs subset of program; calculates BDamage values and writes
                # them to a DataFrame
                pdb.rabdam_dataframe(run='rabdam_dataframe')
            elif vars(args)['output'].lower() in ['analysis']:
                # Runs subset of program; generates output analysis files from
                # pre-calculated BDamage values
                pdb.rabdam_analysis(run='rabdam_analysis')


# Prints total program run time to screen
runtime = time.time() - start
mins = math.floor(runtime/60)
secs = math.fmod(runtime, 60)
if mins == 0:
    print 'Program run time: %02.3f sec\n\n' % secs
elif mins >= 1:
    print 'Program run time: %01.0f min %02.3f sec\n\n' % (mins, secs)
