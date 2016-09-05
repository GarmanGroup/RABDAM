

import numpy as np

def get_xyz_from_objects(auAtomList, atomList):
    # Returns numpy arrays of the x, y and z coordinates of the atoms.

    au_atom_coords = np.zeros([len(auAtomList), 3])
    surrounding_atom_coords = np.zeros([len(atomList), 3])
    for i, atom in enumerate(auAtomList):
        au_atom_coords[i, :] = np.array([atom.xyzCoords[0][0],
                                         atom.xyzCoords[1][0],
                                         atom.xyzCoords[2][0]])
    for j, surr_atom in enumerate(atomList):
        surrounding_atom_coords[j, :] = np.array([surr_atom.xyzCoords[0][0],
                                                  surr_atom.xyzCoords[1][0],
                                                  surr_atom.xyzCoords[2][0]])
    return au_atom_coords, surrounding_atom_coords


def calc_packing_density(xyz_au_atom, xyz_surr_atom, pack_dens_thresh):
    # Calculates the packing density of each atom in the subset of atoms to
    # be considered for B_damage analysis.

    num_au_atoms = xyz_au_atom.shape[0]
    packing_density_array = np.zeros([num_au_atoms])

    for i in xrange(0, num_au_atoms):
        distances = np.sqrt(np.square(xyz_surr_atom - xyz_au_atom[i, :]).sum(axis=1))
        packing_density_array[i] = np.sum(distances < pack_dens_thresh) - 1  # Subtract 1 to exclude the atom counting itself.

    return packing_density_array


def write_pckg_dens_to_atoms(au_atoms, packing_density_array):
    # Writes packing density values to their corresponding atom objects.

    for i, atom in enumerate(au_atoms):
        atom.pd = int(packing_density_array[i])


def calcBdam(bdamatomList, window):
    # All atoms to be considered for B_damage analysis are ordered via their
    # packing density values; the B_damage value of each atom is then
    # calculated as the ratio of its B factor as compared to the average of
    # the B factor values of similarly (identified via sliding window)
    # packed atoms.

    import pandas as pd
    import math

    # Initialises lists of the atom properties required to calculate B_damage.
    ATMNUM = [None]*len(bdamatomList)
    BFAC = [None]*len(bdamatomList)
    PD = [None]*len(bdamatomList)

    # Lists are filled with property values associated with each of the atoms
    # considered for B_damage analysis, then concatenated into the columns
    # of a DataFrame.
    for index, atm in enumerate(bdamatomList):
        ATMNUM[index] = atm.atomNum
        BFAC[index] = atm.bFactor
        PD[index] = atm.pd

    df = pd.DataFrame({'ATMNUM': ATMNUM,
                       'BFAC': BFAC,
                       'PD': PD})

    # DataFrame rows are sorted by packing density (and next by atom number
    # in cases of equal packing density). The B factor values are then used
    # to calculate the average B factor values (calculated via a rolling mean
    # approach using a window size as specified in INPUT.txt (default = 2%))
    # for use in B_damage calculations. In the cases of those atoms which lie
    # too close to either edge of the packing density distribution to lie at
    # the centre of a full-sized window, the average B factor value of each
    # of these atoms is taken from the closest complete window.
    df = df.sort_values(by=['PD', 'ATMNUM'], ascending=[True, True])
    df = df.reset_index(drop=True)

    ser = df[df.columns[1]]
    ser = ser.rename('AVRG_BF')
    ser = ser.rolling(window=window, center=True).mean()
    ser = ser.fillna(0)

    index_list = range(0, len(bdamatomList))
    index = pd.Series(index_list)
    index = index.rename('INDEX')

    df = pd.concat([df, ser, index], axis=1)

    df.loc[(df['AVRG_BF'] == 0) & (df['INDEX'] <= (math.floor(window)-1)),
           'AVRG_BF'] += df.BFAC.values[0:window].mean(axis=0)
    df.loc[(df['AVRG_BF'] == 0) & (df['INDEX'] >= (len(bdamatomList) - math.floor(window / 2))),
           'AVRG_BF'] += df.BFAC.values[(len(bdamatomList) - window):len(bdamatomList)].mean(axis=0)

    df = df.sort_values(by='ATMNUM', ascending=True)
    df = df.reset_index(drop=True)

    # The B_damage value of each atom in the DataFrame is calculated as the
    # ratio of its B factor value to its associated average B factor value.
    for index, value in enumerate(df.ATMNUM.values):
        for atm in bdamatomList:
            if value == atm.atomNum:
                atm.avrg_bf = float(df.AVRG_BF.values[index])
                atm.bd = atm.bFactor / float(df.AVRG_BF.values[index])
