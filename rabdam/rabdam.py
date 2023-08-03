
# RABDAM
# Copyright (C) 2023 Garman Group, University of Oxford

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

class ArgumentError(Exception):
    pass


class FileDoesNotExistError(Exception):
    pass


def parse_command_line_arguments(command_line, test=False):
    """
    Reads in command line arguments and checks that there are no unexpected
    argument values.
    """

    import argparse
    import os
    import sys

    if test is True:
        if sys.version_info[0] < 3:
            from Subroutines.checkDependencies import check_RABDAM_dependencies
        else:
            from rabdam.Subroutines.checkDependencies import check_RABDAM_dependencies
    else:
        if __name__ == '__main__':
            from Subroutines.checkDependencies import check_RABDAM_dependencies
        else:
            if sys.version_info[0] < 3:
                from Subroutines.checkDependencies import check_RABDAM_dependencies
            else:
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
    args = parser.parse_args(command_line)

    # Checks system for RABDAM dependencies
    if vars(args)['dependencies']:
        check_RABDAM_dependencies()
        sys.exit()

    # If specified, checks that file path to input file exists
    if vars(args)['input'] is not None:
        input_file_path = vars(args)['input']
        if not os.path.isfile(input_file_path):
            raise FileDoesNotExistError(
                'Specified input file {} does not exist'.format(input_file_path)
            )

    # Checks options specified for output files are allowed. If no output files
    # are specified, sets default option (= generate all output files).
    output_options = ['csv', 'pdb', 'cif', 'kde', 'bnet', 'summary']
    if vars(args)['output'] is None:
        vars(args)['output'] = output_options
    vars(args)['output'] = [item.lower() for item in vars(args)['output']]

    unrecognised_options = []
    for option in vars(args)['output']:
        if not option in output_options:
            unrecognised_options.append(option)
    if len(unrecognised_options) >= 1:
        raise ArgumentError(
            'Unrecognised output options {} requested.\nPlease specify one or more '
            'of the following: {}.\n'.format(unrecognised_options, output_options)
        )

    # Checks options specified for RABDAM run are allowed. If not specified, by
    # default runs the full analysis pipeline.
    all_run_options = ['dataframe', 'df', 'analysis', None]
    if not vars(args)['run'] in all_run_options:
        raise ArgumentError(
            'Unrecognised run option: {}\nPlease either specify one of the '
            'following options: {}\nOR remove the -r flag (which will direct '
            'RABDAM to run the full analysis pipeline).'.format(vars(args)['run'],
             all_run_options)
        )
    if not vars(args)['run'] is None:
        vars(args)['run'] = vars(args)['run'].lower()

    return args


def parse_input_file_arguments(splitArgs):
    """
    Reads in program parameters listed in user-specified input file
    """

    import os

    # Initialises default program options
    cwd = os.getcwd()
    outputLoc = cwd
    batchVal = False
    overwriteVal = False
    pdtVal = 7.0
    windowVal = 0.02
    protOrNAVal = 'protein'
    hetatmVal = False
    removeAtomsList = []
    addAtomsList = []
    highlightAtomsList = []
    orig_pdb = False
    au_pdb = False
    uc_pdb = False
    auc_pdb = False
    ta_pdb = False

    # If input file is provided, program options are updated to the values
    # the user has specified in it
    for x in range(0, len(splitArgs)):
        # Specifies location to which output 'Logfiles' directory is written
        if splitArgs[x][0:9].lower() == 'outputdir':
            outputLoc = splitArgs[x].split('=')[-1]
            if outputLoc == '':
                outputLoc = cwd
            outputLoc = outputLoc.replace('\\', '/')
            if not os.path.isdir(outputLoc):
                raise FileDoesNotExistError(
                    'Directory {} does not exist'.format(outputLoc)
                )

        # Specifies when an error is encountered whether to exit the program
        # (default) or to continue to the next listed structure
        elif splitArgs[x][0:13].lower() == 'batchcontinue':
            batchVal = splitArgs[x].split('=')[-1].lower()
            if batchVal in ['true', 'yes', 't', 'y']:
                batchVal = True
            elif batchVal in ['false', 'no', 'f', 'n']:
                batchVal = False
            else:
                raise ArgumentError(
                    'Unrecognised value for batchcontinue: {}'.format(batchVal)
                )

        # Specifies whether or not to overwrite files of the same name as the
        # new output files to be created (default = do not overwrite)
        elif splitArgs[x][0:9].lower() == 'overwrite':
            overwriteVal = splitArgs[x].split('=')[-1].lower()
            if overwriteVal in ['true', 'yes', 't', 'y']:
                overwriteVal = True
            elif overwriteVal in ['false', 'no', 'f', 'n']:
                overwriteVal = False
            else:
                raise ArgumentError(
                    'Unrecognised value for overwrite: {}'.format(overwriteVal)
                )

            # If command line output is directed to a file, sets overwriteVal
            # to True to prevent program from requiring user input to decide
            # whether to overwrite any pre-existing file(s) with the same file
            # path(s) as the output file(s) to be written
            if os.fstat(0) != os.fstat(1):
                overwriteVal = True

        # Specifies packing density threshold (default = 7)
        elif splitArgs[x][0:3].lower() == 'pdt':
            pdtVal = splitArgs[x].split('=')[-1]
            try:
                if float(pdtVal) != int(pdtVal):
                    raise ArgumentError(
                        'Value provided for packing density is not an integer: '
                        '{}'.format(pdtVal)
                    )
                else:
                    pdtVal = float(pdtVal)
            except ValueError:
                raise ArgumentError(
                    'Value provided for packing density is not an integer: '
                    '{}'.format(pdtVal)
                )

        # Specifies size of sliding window (as a percentage of the total number
        # of atoms considered for BDamage analysis)
        elif splitArgs[x][0:10].lower() == 'windowsize':
            windowVal = splitArgs[x].split('=')[-1]
            try:
                windowVal = float(windowVal)
            except ValueError:
                raise ArgumentError(
                    'Value provided for size of sliding window is not a number:'
                    ' {}'.format(windowVal)
                )
            if not 0 < windowVal < 1:
                raise ArgumentError(
                    'Value provided for size of sliding window is not in the '
                    'range 0 - 1: {}'.format(windowVal)
                )

        # Specifies whether to include protein atoms alone, nucleic atoms
        # alone, or both atom types (if present) in the BDamage calculation
        elif splitArgs[x][0:20].lower() == 'proteinornucleicacid':
            protOrNAVal = splitArgs[x].split('=')[-1].lower()
            if not protOrNAVal in ['protein', 'na', 'nucleicacid', 'proteinna']:
                raise ArgumentError(
                    'Unrecognised value for proteinornucleicacid: '
                    '{}'.format(protOrNAVal)
                )

        # Specifies whether to remove HETATM from the BDamage calculation or to
        # retain them
        elif splitArgs[x][0:6].lower() == 'hetatm':
            hetatmVal = splitArgs[x].split('=')[-1].lower()
            if hetatmVal == 'keep':
                hetatmVal = True
            elif hetatmVal == 'remove':
                hetatmVal = False
            else:
                raise ArgumentError(
                    'Unrecognised value for hetatm: {}'.format(hetatmVal)
                )

        # Lists atoms (via either their atom numbers or their residue names) to
        # be removed from the BDamage calculation. (This is useful to remove
        # additional atoms not covered by the proteinOrNucleicAcid / HETATM
        # options.)
        elif splitArgs[x][0:11].lower() == 'removeatoms':
            removeAtomsArg = splitArgs[x].split('=')[-1].upper()
            removeAtomsList = []
            if removeAtomsArg != '':
                removeAtomsArg = removeAtomsArg.split(';')
                for item in removeAtomsArg:
                    removeAtomsRange = item.split('-')
                    if len(removeAtomsRange) == 2:
                        try:
                            min_val = int(removeAtomsRange[0])
                            max_val = int(removeAtomsRange[1])
                            removeAtomsRange = range(min_val, (max_val+1))
                            for number in removeAtomsRange:
                                removeAtomsList.append(str(number))
                        except ValueError:
                            raise ArgumentError(
                                'Unrecognised input: {} for remove atoms - if '
                                'input contains "-", expecting a numeric '
                                'range'.format(removeAtomsRange)
                            )
                    elif len(removeAtomsRange) == 1:
                        removeAtomsList.append(str(removeAtomsRange[0]))
                    else:
                        raise ArgumentError(
                            'Unrecognised input: {} for remove atoms - if '
                            'input contains "-", expecting a numeric '
                            'range'.format(removeAtomsRange)
                        )

        # Lists atoms (via either their atom numbers or their residue names) to
        # be included in the BDamage calculation. (This is useful to add back
        # in a subset of atoms of interest that has been removed by the
        # proteinOrNucleicAcid / HETATM / removeAtoms options.)
        elif splitArgs[x][0:8].lower() == 'addatoms':
            addAtomsArg = splitArgs[x].split('=')[-1].upper()
            addAtomsList = []
            if addAtomsArg != '':
                addAtomsArg = addAtomsArg.split(';')
                for item in addAtomsArg:
                    addAtomsRange = item.split('-')
                    if len(addAtomsRange) == 2:
                        try:
                            min_val = int(addAtomsRange[0])
                            max_val = int(addAtomsRange[1])
                            addAtomsRange = range(min_val, (max_val+1))
                            for number in addAtomsRange:
                                addAtomsList.append(str(number))
                        except ValueError:
                            raise ArgumentError(
                                'Unrecognised input: {} for add atoms - if '
                                'input contains "-", expecting a numeric '
                                'range'.format(addAtomsRange)
                            )
                    elif len(addAtomsRange) == 1:
                        addAtomsList.append(str(addAtomsRange[0]))
                    else:
                        raise ArgumentError(
                            'Unrecognised input: {} for add atoms - if '
                            'input contains "-", expecting a numeric '
                            'range'.format(addAtomsRange)
                        )

        # Lists atoms (via their atom numbers) whose BDamage values are to be
        # indicated on the kernel density estimate of the BDamage distribution
        # output by the program. Note that it is recommended no more than 6
        # atoms are listed in the highlightAtoms option in the input file
        # (beyond 6 atoms, the colour scheme will repeat itself, and in
        # addition the key may not fit onto the graph).
        elif splitArgs[x][0:14].lower() == 'highlightatoms':
            highlightAtomsArg = splitArgs[x].split('=')[-1].upper()
            highlightAtomsList = []
            if highlightAtomsArg != '':
                highlightAtomsArg = highlightAtomsArg.split(';')
                for item in highlightAtomsArg:
                    highlightAtomsRange = item.split('-')
                    if len(highlightAtomsRange) == 2:
                        try:
                            min_val = int(highlightAtomsRange[0])
                            max_val = int(highlightAtomsRange[1])
                            highlightAtomsRange = range(min_val, (max_val+1))
                            for number in highlightAtomsRange:
                                highlightAtomsList.append(str(number))
                        except ValueError:
                            raise ArgumentError(
                                'Unrecognised input: {} for highlight atoms - '
                                'if input contains "-", expecting a numeric '
                                'range'.format(highlightAtomsRange)
                            )
                    elif len(highlightAtomsRange) == 1:
                        highlightAtomsList.append(str(highlightAtomsRange[0]))
                    else:
                        raise ArgumentError(
                            'Unrecognised input: {} for highlight atoms - if '
                            'input contains "-", expecting a numeric '
                            'range'.format(highlightAtomsRange)
                        )

        # Specifies whether to save the input PDB file fed into the program
        elif splitArgs[x][0:13].lower() == 'createorigpdb':
            orig_pdb = splitArgs[x].split('=')[-1].lower()
            if orig_pdb in ['true', 'yes', 't', 'y']:
                orig_pdb = True
            elif orig_pdb in ['false', 'no', 'f', 'n']:
                orig_pdb = False
            else:
                raise ArgumentError(
                    'Unrecognised value for createorigpdb: {}'.format(orig_pdb)
                )

        # Specifies whether to create a PDB file of the (filtered) asymmetric
        # unit
        elif splitArgs[x][0:11].lower() == 'createaupdb':
            au_pdb = splitArgs[x].split('=')[-1].lower()
            if au_pdb in ['true', 'yes', 't', 'y']:
                au_pdb = True
            elif au_pdb in ['false', 'no', 'f', 'n']:
                au_pdb = False
            else:
                raise ArgumentError(
                    'Unrecognised value for createaupdb: {}'.format(au_pdb)
                )

        # Specifies whether to create a PDB file of the unit cell
        elif splitArgs[x][0:11].lower() == 'createucpdb':
            uc_pdb = splitArgs[x].split('=')[-1].lower()
            if uc_pdb in ['true', 'yes', 't', 'y']:
                uc_pdb = True
            elif uc_pdb in ['false', 'no', 'f', 'n']:
                uc_pdb = False
            else:
                raise ArgumentError(
                    'Unrecognised value for createucpdb: {}'.format(uc_pdb)
                )

        # Specifies whether to create a PDB file of the 3x3x3 unit cell
        # assembly
        elif splitArgs[x][0:12].lower() == 'createaucpdb':
            auc_pdb = splitArgs[x].split('=')[-1].lower()
            if auc_pdb in ['true', 'yes', 't', 'y']:
                auc_pdb = True
            elif auc_pdb in ['false', 'no', 'f', 'n']:
                auc_pdb = False
            else:
                raise ArgumentError(
                    'Unrecognised value for createaucpdb: {}'.format(auc_pdb)
                )

        # Specifies whether to create a PDB file of the trimmed atoms assembly
        # (which consists of the asymmetric unit, plus every atom in the
        # 3x3x3 unit cell assembly within a PDT (default=7) Angstrom radius of
        # the asymmetric unit)
        elif splitArgs[x][0:11].lower() == 'createtapdb':
            ta_pdb = splitArgs[x].split('=')[-1].lower()
            if ta_pdb in ['true', 'yes', 't', 'y']:
                ta_pdb = True
            elif ta_pdb in ['false', 'no', 'f', 'n']:
                ta_pdb = False
            else:
                raise ArgumentError(
                    'Unrecognised value for createtapdb: {}'.format(ta_pdb)
                )

        else:
            if splitArgs[x] != '':
                raise ArgumentError(
                    'Unrecognised argument {}'.format(splitArgs[x])
                )

    input_arguments = {'outputDir': outputLoc,
                       'batchRun': batchVal,
                       'overwrite': overwriteVal,
                       'PDT': pdtVal,
                       'windowSize': windowVal,
                       'protOrNA': protOrNAVal,
                       'HETATM': hetatmVal,
                       'removeAtoms': removeAtomsList,
                       'addAtoms': addAtomsList,
                       'highlightAtoms': highlightAtomsList,
                       'createOrigpdb': orig_pdb,
                       'createAUpdb': au_pdb,
                       'createUCpdb': uc_pdb,
                       'createAUCpdb': auc_pdb,
                       'createTApdb': ta_pdb}

    return input_arguments


def main(test=False):
    """
    Runs RABDAM pipeline
    """

    import math
    import sys
    import time

    if __name__ == '__main__':
        from Subroutines.CalculateBDamage import rabdam
    else:
        if sys.version_info[0] < 3:
            from Subroutines.CalculateBDamage import rabdam
        else:
            from rabdam.Subroutines.CalculateBDamage import rabdam

    # Parses in command line arguments
    args = parse_command_line_arguments(sys.argv[1:])

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

    # Reads in program options from input file specified by -i flag
    if vars(args)['input'] is not None:
        input_file_path = vars(args)['input']
        input_file_path = input_file_path.replace('\\', '/')
        fileCont = open(input_file_path, 'r')
        functionArgs = fileCont.read()
        fileCont.close()
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

    # Parses in input arguments if input file specified with -i flag / sets to
    # default if only file name / accession code specified with -f flag
    input_arguments = parse_input_file_arguments(splitArgs)

    # Runs the BDamage calculation for every specified PDB file
    for item in pathToInputList:
        # Initialises rabdam object
        pdb = rabdam(
            pathToInput=item,
            outputDir=input_arguments['outputDir'],
            batchRun=input_arguments['batchRun'],
            overwrite=input_arguments['overwrite'],
            PDT=input_arguments['PDT'],
            windowSize=input_arguments['windowSize'],
            protOrNA=input_arguments['protOrNA'],
            HETATM=input_arguments['HETATM'],
            removeAtoms=input_arguments['removeAtoms'],
            addAtoms=input_arguments['addAtoms'],
            highlightAtoms=input_arguments['highlightAtoms'],
            createOrigpdb=input_arguments['createOrigpdb'],
            createAUpdb=input_arguments['createAUpdb'],
            createUCpdb=input_arguments['createUCpdb'],
            createAUCpdb=input_arguments['createAUCpdb'],
            createTApdb=input_arguments['createTApdb']
        )

        if vars(args)['run'] in [None, 'dataframe', 'df']:
            # Calculates BDamage values and writes them to a DataFrame
            pdb.rabdam_dataframe(test)
        if vars(args)['run'] in [None, 'analysis']:
            # Generates output analysis files from pre-calculated BDamage values
            pdb.rabdam_analysis(output_options=vars(args)['output'])

    # Prints total program run time to screen
    runtime = time.time() - start
    mins = math.floor(runtime/60)
    secs = math.fmod(runtime, 60)
    if mins == 0:
        print('Program run time: %02.3f sec\n' % secs)
    elif mins >= 1:
        print('Program run time: %01.0f min %02.3f sec\n' % (mins, secs))

    print('Please cite:\nShelley KL, Dixon TPE, Brooks-Bartlett JC & Garman EF'
          ' (2018). J Appl Cryst 51, 552-559.\n\n')
    print('If using BDamage, please cite:\nGerstel M, Deane CM & Garman EF '
          '(2015).  J Synchrotron Radiat 22, 201-212.\n\n')
    print('If using Bnet, please cite:\nShelley KL & Garman EF (2022). '
          'Nat Commun 13, 1314.\n\n')

# Runs 'main' function if rabdam.py is run as a script
if __name__ == '__main__':
    main()
