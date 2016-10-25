

def genPDBCURinputs(PDBCURinputFile, asymmetricUnit):
    # Creates input file for PDBCUR program specifying the operations to be
    # performed, namely:
    # - delete hydrogen atoms
    # - remove all zero occupancy atoms
    # - retain only highest occupancy conformations (if equal occupancies then
    #   only one (the first listed) is selected), and set occupancy values of
    #   retained atoms to 1
    # - remove anisotropic B factors
    # - if specified create a unit cell from the asymmetric unit plus its
    #   associated symmetry operations

    print 'Creating input file for PDBCUR at %s' % PDBCURinputFile
    with open(PDBCURinputFile, 'w') as f:
        f.write('delhydrogen\n')
        f.write('cutocc\n')
        f.write('mostprob\n')
        f.write('noanisou\n')
        if asymmetricUnit is False:
            f.write('genunit\n')
        f.close


def runPDBCUR(pathToPDB, PDBCURoutputPDB, PDBCURinputFile, PDBCURlog):
    # Runs PDBCUR from the command line with the operations specified in
    # PDBCUR input file.

    import os

    runPDBCURcommand = 'pdbcur xyzin %s xyzout %s < %s > %s' % (
        pathToPDB, PDBCURoutputPDB, PDBCURinputFile, PDBCURlog
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
