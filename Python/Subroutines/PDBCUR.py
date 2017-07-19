
def clean_pdb_file(pathToPDB, pdb_file_path):
    import math
    import pandas as pd
    from parsePDB import atom

    orig_pdb_lines = []
    header_lines = []
    footer_lines = []
    unit_cell_params = []
    orig_pdb = open('%s' % pathToPDB, 'r')
    header = True
    footer = False
    for line in orig_pdb:
        if line[0:6].strip() == 'CRYST1':
            params = line.split()
            a = float(params[1])
            b = float(params[2])
            c = float(params[3])
            alpha = math.radians(float(params[4]))
            beta = math.radians(float(params[5]))
            gamma = math.radians(float(params[6]))
            unit_cell_params.extend((a, b, c, alpha, beta, gamma))
        elif line[0:6].strip() in ['ATOM', 'HETATM']:  # Removes anisotropic B factors
            orig_pdb_lines.append(line)
            header = False
        elif line[0:6].strip() == 'CONECT':
            footer = True

        if header is True:
            header_lines.append(line)
        elif footer is True:
            footer_lines.append(line)
    orig_pdb.close()

    # Removes hydrogens and 0 occupancy atoms
    filtered_pdb_lines = [line for line in orig_pdb_lines if line[76:78].strip()
                          != 'H' and float(line[54:60].strip()) != 0.0]

    # Removes alternate conformers
    alternate_conformers_reschainnum = []
    alternate_conformers_label = []
    alternate_conformers_occupancy = []
    for line in filtered_pdb_lines:
        if line[16:17].strip() != '':
            alternate_conformers_reschainnum.append(line[21:26].replace(' ', ''))
            alternate_conformers_label.append(line[16:17].strip())
            alternate_conformers_occupancy.append(line[54:60].strip())

    df = pd.DataFrame({'reschainnum': alternate_conformers_reschainnum,
                       'conformer': alternate_conformers_label,
                       'occupancy': alternate_conformers_occupancy})
    df = df.drop_duplicates()
    reschainnum = df['reschainnum'].tolist()
    conformer = df['conformer'].tolist()
    occupancy = df['occupancy'].tolist()

    alternate_conformers = {}
    reschainnum_set = set(reschainnum)
    for number in reschainnum_set:
        indices = []
        a = 'A1'
        b = 'A'
        c = 0
        for index_1, value_1 in enumerate(reschainnum):
            if value_1 == number:
                indices.append(index_1)
        for index_2 in indices:
            if occupancy[index_2] > c:
                a = reschainnum[index_2]
                b = conformer[index_2]
                c = occupancy[index_2]
        alternate_conformers[a] = b

    clean_au_list = []
    clean_au_file = open('%s_asymmetric_unit.pdb' % pdb_file_path, 'w')
    for line in header_lines:
        clean_au_file.write(line)
    for line in filtered_pdb_lines:
        if line[16:17].strip() == '':
            clean_au_file.write(line)
            y = atom()
            y.lineID = str(line[0:6].strip())
            y.atomNum = int(line[6:11].strip())
            y.atomType = str(line[12:16].strip())
            y.resiType = str(line[17:20].strip())
            y.chainID = str(line[21:22].strip())
            y.resiNum = int(line[22:26].strip())
            y.xyzCoords = [[float(line[30:38].strip())],
                           [float(line[38:46].strip())],
                           [float(line[46:54].strip())]]
            y.occupancy = float(line[54:60].strip())
            y.bFactor = float(line[60:66].strip())
            y.atomID = str(line[76:78].strip())
            y.charge = str(line[78:80].strip())
            clean_au_list.append(y)
        elif line[16:17].strip() != '':
            if line[21:26].replace(' ', '') in alternate_conformers:
                if alternate_conformers[line[21:26].replace(' ', '')] == line[16:17].strip():
                    clean_au_file.write(line)
                    y = atom()
                    y.lineID = str(line[0:6].strip())
                    y.atomNum = int(line[6:11].strip())
                    y.atomType = str(line[12:16].strip())
                    y.resiType = str(line[17:20].strip())
                    y.chainID = str(line[21:22].strip())
                    y.resiNum = int(line[22:26].strip())
                    y.xyzCoords = [[float(line[30:38].strip())],
                                   [float(line[38:46].strip())],
                                   [float(line[46:54].strip())]]
                    y.occupancy = float(line[54:60].strip())
                    y.bFactor = float(line[60:66].strip())
                    y.atomID = str(line[76:78].strip())
                    y.charge = str(line[78:80].strip())
                    clean_au_list.append(y)
    for line in footer_lines:
        clean_au_file.write(line)
    clean_au_file.close()

    clean_au_file_name = '%s_asymmetric_unit.pdb' % pdb_file_path

    return (clean_au_file_name, clean_au_list, header_lines, footer_lines,
            unit_cell_params)


def genPDBCURinputs(PDBCURinputFile):
    # Creates input file for the CCP4 suite program PDBCUR instructing it to
    # create the unit cell from an input pdb file of an asymmetric unit

    print 'Creating input file for PDBCUR at %s' % PDBCURinputFile
    input_file = open(PDBCURinputFile, 'w')
    input_file.write('genunit\n')
    input_file.close()


def runPDBCUR(clean_au_file, PDBCURoutputPDB, PDBCURinputFile, PDBCURlog):
    # Runs PDBCUR from the command line with the operations specified in
    # PDBCUR input file.

    import os

    runPDBCURcommand = 'pdbcur xyzin %s xyzout %s < %s > %s' % (
        clean_au_file, PDBCURoutputPDB, PDBCURinputFile, PDBCURlog
        )
    print 'Running PDBCUR (Winn et al. 2011) to process the PDB file'
    os.system(runPDBCURcommand)
    print 'PDBCUR log is printed below\n'
    PDBCURlogText = open(PDBCURlog, 'r')
    for line in PDBCURlogText:
        print line
    PDBCURlogText.close()
    os.remove(PDBCURlog)
    os.remove(PDBCURinputFile)

    print('\n################################################################\n'
          '\n################################################################\n'
          '\n################################################################\n')
