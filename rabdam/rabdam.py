
# RABDAM
# Copyright (C) 2018 Garman Group, University of Oxford

# This file is part of RABDAM.

# RABDAM is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation, either version 3 of
# the License, or (at your option) any later version.

# RABDAM is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General
# Public License along with this program.  If not, see
# <http://www.gnu.org/licenses/>.

# An outer layer to the pipeline scripts. Depending upon the flags specified
# in the command line input, this script will run either the complete / a
# subsection of the pipeline.

def main():
    import sys
    import os
    import argparse
    import time
    import math
    import numpy as np
    import argparse

    if __name__ == '__main__' or 'CCP4' in list(os.environ.keys()):
        from Subroutines.CalculateBDamage import rabdam
        from Subroutines.checkDependencies import check_RABDAM_dependencies
    else:
        from rabdam.Subroutines.CalculateBDamage import rabdam
        from rabdam.Subroutines.checkDependencies import check_RABDAM_dependencies

    # Reads in command line inputs. There are three recognised flags: -i, -f
    # and -o. Program inputs are specified by either the -i or the -f flag;
    # provision of one of these flags is compulsory. Program outputs are
    # specified by the optional -r and -o flags.
    parser = argparse.ArgumentParser()
    input_file_group = parser.add_mutually_exclusive_group(required=True)
    input_file_group.add_argument('--dependencies', action='store_true',
                                  help='Checks whether your system has the '
                                  'necessary packages / programs installed in '
                                  'order to be able to run RABDAM')
    input_file_group.add_argument('-i', '--input', help='Absolute path to '
                                  'input file listing program parameter values')
    input_file_group.add_argument('-f', '--pdb_or_mmcif_file', nargs='+',
                                  help='Specifies input pdb file (via either '
                                  'a 4 character PDB accession code or an '
                                  'absolute file path) for BDamage analysis - '
                                  'this option allows the RABDAM program to '
                                  'be run (using default program parameter '
                                  'values) without providing an input file '
                                  'listing program options')
    parser.add_argument('-r', '--run', help='Specifies whether to run the '
                        'complete program (= default), to calculate BDamage '
                        'values only ("df" / "dataframe"), or to analyse '
                        'pre-calculated BDamage values only ("analysis")')
    parser.add_argument('-o', '--output', nargs='+', help='Specifies the '
                        'output files to write (default = all files written)')
    args = parser.parse_args()

    # Checks system for RABDAM dependencies
    if vars(args)['dependencies']:
        check_RABDAM_dependencies()
        sys.exit()

    # Initialises program start time
    start = time.time()
    startIndex = time.localtime()
    year = startIndex.tm_year
    month = startIndex.tm_mon
    day = startIndex.tm_mday
    hour = startIndex.tm_hour
    minute = startIndex.tm_min
    second = startIndex.tm_sec

    print('\nThis program was run on %d/%d/%d at %02.0f:%02.0f:%02.0f\n\n' % (
        day, month, year, hour, minute, second
        ))

    # Sets default option for -o flag (default = generate all output files bar
    # the summary file)
    if vars(args)['output'] is None:
        vars(args)['output'] = ['csv', 'bdam', 'kde', 'bnet', 'summary']
    output_options = [item.lower() for item in vars(args)['output']]

    cwd = os.getcwd()
    # Reads in program options from input file specified by -i flag
    if vars(args)['input'] is not None:
        # Reads in PDB and cif file name(s) listed in input file
        input_file_loc = vars(args)['input']
        input_file_loc = input_file_loc.replace('\\', '/')
        input_file_loc_list = input_file_loc.split('/')
        input_file = input_file_loc_list[len(input_file_loc_list)-1]
        if len(input_file_loc_list) > 1:
            input_file_loc = input_file_loc.replace(input_file, '')
            os.chdir(input_file_loc)
        fileCont = open(input_file, 'r')
        functionArgs = fileCont.read()
        fileCont.close()
        os.chdir(cwd)
        functionArgs = functionArgs.replace('\n', '')
        functionArgs = functionArgs.replace('\r', '')
        splitArgs = functionArgs.split(',')
        pathToInputList = [item.strip() for item in splitArgs if '=' not in item]
        # Reads in remaining program options specified in input file
        functionArgs = functionArgs.replace(' ', '')
        splitArgs = functionArgs.split(',')
        splitArgs = [item.strip() for item in splitArgs if '=' in item]

    # Reads in PDB file name(s) specified by -f flag listed in command line input
    elif vars(args)['pdb_or_mmcif_file'] is not None:
        pathToInputList = [item.strip() for item in vars(args)['pdb_or_mmcif_file']]
        splitArgs = []

    # Checks that the user has specified at least one input structure for RABDAM
    # analysis
    if len(pathToInputList) == 0:
        sys.exit('Input file / accession code not specified')

    # Initialises default program options
    outputLoc = cwd
    batchVal = False
    overwriteVal = False
    pdtList = [7.0]
    windowList = [0.02]
    protOrNAVal = 'protein'
    hetatmVal = False
    removeAtomsList = []
    addAtomsList = []
    highlightAtomsList = []
    origVal = False
    auVal = False
    ucVal = False
    aucVal = False
    taVal = False

    # If input file is provided, program options are updated to the values
    # the user has specified in it
    for x in range(0, len(splitArgs)):
        # Specifies location to which output 'Logfiles' directory is written
        if splitArgs[x][0:3].lower() == 'dir':
            dirArg = splitArgs[x].split('=')
            outputLoc = dirArg[len(dirArg)-1]
            if outputLoc == '':
                outputLoc = cwd
            outputLoc = outputLoc.replace('\\', '/')

        # Specifies if an error is encountered when analysing the current
        # structure whether to exit the program (default) or to continue to
        # analyse the next listed structure
        elif splitArgs[x][0:13].lower() == 'batchcontinue':
            batchArg = splitArgs[x].split('=')
            batchVal = batchArg[len(batchArg)-1].lower()
            if batchVal in ['true', 'yes', 't', 'y']:
                batchVal = True
            elif batchVal in ['false', 'no', 'f', 'n']:
                batchVal = False

        # Specifies whether or not to overwrite files of the same name as the
        # new output files to be created
        elif splitArgs[x][0:9].lower() == 'overwrite':
            overwriteArg = splitArgs[x].split('=')
            overwriteVal = overwriteArg[len(overwriteArg)-1].lower()
            if overwriteVal in ['true', 'yes', 't', 'y']:
                overwriteVal = True
            elif overwriteVal in ['false', 'no', 'f', 'n']:
                overwriteVal = False

            # If command line output is directed to a file, sets overwriteVal
            # to True to prevent program from requiring user input to decide
            # whether to overwrite any pre-existing file(s) with the same file
            # path(s) as the output file(s) to be written
            if os.fstat(0) != os.fstat(1):
                overwriteVal = True

        # Specifies packing density threshold
        elif splitArgs[x][0:3].lower() == 'pdt':
            pdtArg = splitArgs[x].split('=')
            pdtArg = pdtArg[len(pdtArg)-1]
            if '-' in pdtArg:
                pdtRange = pdtArg.split('-')
                min_pdt = float(pdtRange[0])
                max_pdt = float(pdtRange[len(pdtRange)-1])
                pdtList = np.arange(min_pdt, max_pdt, 1.0)
            else:
                pdtList = [float(pdtArg)]

        # Specifies size of sliding window (as a percentage of the total number
        # of atoms considered for BDamage analysis)
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

        # Specifies whether to include protein atoms alone, nucleic atoms
        # alone, or both atom types (if present) in the BDamage calculation
        elif splitArgs[x][0:20].lower() == 'proteinornucleicacid':
            protOrNAArg = splitArgs[x].split('=')
            protOrNAVal = protOrNAArg[len(protOrNAArg)-1].lower()
            # Temporary variable reassignment to prevent nucleic acid analysis
            # with RABDAM before the complete functionality has been introduced
            if protOrNAVal in ['na', 'nucleicacid']:
                sys.exit('RABDAM is currently only suitable for assessing '
                         'radiation damage\n'
                         'to the protein component of macromolecular '
                         'structures.\n'
                         'We hope to extend the program to incorporate '
                         'nucleic acid analysis\n'
                         'shortly - in the meantime, please restrict your '
                         'RABDAM analysis to\n'
                         'protein atoms only.')

        # Specifies whether to remove HETATM from the BDamage calculation or to
        # retain them
        elif splitArgs[x][0:6].lower() == 'hetatm':
            hetatmArg = splitArgs[x].split('=')
            hetatmVal = hetatmArg[len(hetatmArg)-1].lower()
            if hetatmVal == 'keep':
                hetatmVal = True
            elif hetatmVal == 'remove':
                hetatmVal = False

        # Lists atoms (via either their atom numbers or their residue names) to
        # be removed from the BDamage calculation. (This is useful to remove
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
                        min_val = int(removeAtomsRange[0])
                        max_val = int(removeAtomsRange[1])
                        removeAtomsRange = range(min_val, (max_val+1))
                        for number in removeAtomsRange:
                            removeAtomsList.append(str(number))
                    elif len(removeAtomsRange) == 1:
                        removeAtomsList.append(removeAtomsRange[0])

        # Lists atoms (via either their atom numbers or their residue names) to
        # be included in the BDamage calculation. (This is useful to add back
        # in a subset of atoms of interest that has been removed by the
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
                        min_val = int(addAtomsRange[0])
                        max_val = int(addAtomsRange[1])
                        addAtomsRange = range(min_val, (max_val+1))
                        for number in addAtomsRange:
                            addAtomsList.append(str(number))
                    elif len(addAtomsRange) == 1:
                        addAtomsList.append(addAtomsRange[0])

        # Lists atoms (via their atom numbers) whose BDamage values are to be
        # indicated on the kernel density estimate of the BDamage distribution
        # output by the program. Note that it is recommended no more than 6
        # atoms are listed in the highlightAtoms option in the input file
        # (beyond 6 atoms, the colour scheme will repeat itself, and in
        # addition the key may not fit onto the graph).
        elif splitArgs[x][0:14].lower() == 'highlightatoms':
            highlightAtomsList = []
            highlightAtomsArg = splitArgs[x].split('=')
            highlightAtomsStr = highlightAtomsArg[len(highlightAtomsArg)-1]
            if highlightAtomsStr != '':
                highlightAtomsSubList = highlightAtomsStr.split(';')
                for item in highlightAtomsSubList:
                    highlightAtomsRange = item.split('-')
                    if len(highlightAtomsRange) == 2:
                        min_val = int(highlightAtomsRange[0])
                        max_val = int(highlightAtomsRange[1])
                        highlightAtomsRange = range(min_val, (max_val+1))
                        for number in highlightAtomsRange:
                            highlightAtomsList.append(str(number))
                    elif len(highlightAtomsRange) == 1:
                        highlightAtomsList.append(highlightAtomsRange[0])

        # Specifies whether to save the input PDB file fed into the program
        elif splitArgs[x][0:11].lower() == 'createorigpdb':
            origArg = splitArgs[x].split('=')
            origVal = origArg[len(origArg)-1].lower()
            if origVal in ['true', 'yes', 't', 'y']:
                origVal = True
            elif origVal in ['false', 'no', 'f', 'n']:
                origVal = False

        # Specifies whether to create a PDB file of the (filtered) asymmetric
        # unit
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

        # Specifies whether to create a PDB file of the 3x3x3 unit cell
        # assembly
        elif splitArgs[x][0:12].lower() == 'createaucpdb':
            aucArg = splitArgs[x].split('=')
            aucVal = aucArg[len(aucArg)-1].lower()
            if aucVal in ['true', 'yes', 't', 'y']:
                aucVal = True
            elif aucVal in ['false', 'no', 'f', 'n']:
                aucVal = False

        # Specifies whether to create a PDB file of the trimmed atoms assembly
        # (which consists of the asymmetric unit, plus every atom in the
        # 3x3x3 unit cell assembly within a PDT (default=7) Angstrom radius of
        # the asymmetric unit)
        elif splitArgs[x][0:11].lower() == 'createtapdb':
            taArg = splitArgs[x].split('=')
            taVal = taArg[len(taArg)-1].lower()
            if taVal in ['true', 'yes', 't', 'y']:
                taVal = True
            elif taVal in ['false', 'no', 'f', 'n']:
                taVal = False

    # Runs the BDamage calculation for every specified PDB file
    for item in pathToInputList:
        for windowVal in windowList:
            for pdtVal in pdtList:
                # Initialises rabdam object
                pdb = rabdam(pathToInput=item, outputDir=outputLoc,
                             batchRun=batchVal, overwrite=overwriteVal,
                             PDT=pdtVal, windowSize=windowVal,
                             protOrNA=protOrNAVal, HETATM=hetatmVal,
                             addAtoms=addAtomsList,
                             removeAtoms=removeAtomsList,
                             highlightAtoms=highlightAtomsList,
                             createOrigpdb=origVal, createAUpdb=auVal,
                             createUCpdb=ucVal, createAUCpdb=aucVal,
                             createTApdb=taVal)
                if vars(args)['run'] is None:
                    # Runs full program
                    pdb.rabdam_dataframe(run='rabdam')
                    pdb.rabdam_analysis(run='rabdam',
                                        output_options=output_options)
                elif vars(args)['run'].lower() in ['dataframe', 'df']:
                    # Runs subset of program; calculates BDamage values and
                    # writes them to a DataFrame
                    pdb.rabdam_dataframe(run='rabdam_dataframe')
                elif vars(args)['run'].lower() in ['analysis']:
                    # Runs subset of program; generates output analysis files
                    # from pre-calculated BDamage values
                    pdb.rabdam_analysis(run='rabdam_analysis',
                                        output_options=output_options)


    # Prints total program run time to screen
    runtime = time.time() - start
    mins = math.floor(runtime/60)
    secs = math.fmod(runtime, 60)
    if mins == 0:
        print('Program run time: %02.3f sec\n' % secs)
    elif mins >= 1:
        print('Program run time: %01.0f min %02.3f sec\n' % (mins, secs))

    print('Please cite:\nShelley, K. L., Dixon, T. P. E., Brooks-Bartlett, '
          'J. C. & Garman, E. F. (2018). J. Appl. Cryst. 51, 552–559\n\n')

# Runs 'main' function if rabdam.py is run as a script
if __name__ == '__main__':
    main()
