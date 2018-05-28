
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

import numpy as np


def get_xyz_from_objects(bdamAtomList):
    # Returns numpy arrays of the x, y and z coordinates of the atoms to be
    # included in the BDamage calculation.

    au_atom_coords = np.zeros([len(bdamAtomList), 3])
    for i, atom in enumerate(bdamAtomList):
        au_atom_coords[i, :] = np.array([atom.xyzCoords[0][0],
                                         atom.xyzCoords[1][0],
                                         atom.xyzCoords[2][0]])

    return au_atom_coords


def calc_packing_density(xyz_au_atom, xyz_surr_atom, pack_dens_thresh):
    # Calculates the packing density of each atom in the subset of atoms to
    # be considered for BDamage analysis.

    num_au_atoms = xyz_au_atom.shape[0]
    packing_density_array = np.zeros([num_au_atoms, 1])

    for i in range(num_au_atoms):
        distances = np.sqrt(np.square(xyz_surr_atom[:, :] - xyz_au_atom[i, :]).sum(axis=1))
        packing_density_array[i][0] = np.sum(distances < pack_dens_thresh) - 1  # Subtract
        # 1 to correct for the atom itself being counted.

    return packing_density_array


def write_pckg_dens_to_atoms(bdamAtomList, packing_density_array):
    # Writes packing density values to their corresponding atom objects.

    for i, atom in enumerate(bdamAtomList):
        atom.pd = packing_density_array[i][0]


def calcBDam(bdamAtomList, window):
    # All atoms to be considered for BDamage analysis are ordered via their
    # packing density values; the BDamage value of each atom is then
    # calculated as the ratio of its B-factor as compared to the average of the
    # B-factor values of similarly (identified via sliding window) packed atoms.

    import pandas as pd
    import math

    # Initialises lists of the atom properties required to calculate BDamage.
    ATMNUM = [None]*len(bdamAtomList)
    BFAC = [None]*len(bdamAtomList)
    PD = [None]*len(bdamAtomList)

    # Lists are filled with property values associated with each of the atoms
    # considered for BDamage analysis, then concatenated into the columns
    # of a DataFrame.
    for index, atm in enumerate(bdamAtomList):
        ATMNUM[index] = atm.atomNum
        BFAC[index] = atm.bFactor
        PD[index] = atm.pd

    df = pd.DataFrame({'ATMNUM': ATMNUM,
                       'BFAC': BFAC,
                       'PD': PD})

    # DataFrame rows are sorted by packing density (and next by atom number
    # in cases of equal packing density). Average B-factor values are then
    # calculated via a rolling mean approach with a window size as specified in
    # the input file (default = 2%). In the cases of those atoms which lie too
    # close to either edge of the packing density distribution to lie at the
    # centre of a full-sized window, the average B-factor value of each of
    # these atoms is taken from the closest complete window.
    df = df.sort_values(by=['PD', 'ATMNUM'], ascending=[True, True])
    df = df.reset_index(drop=True)

    ser = df[df.columns[1]]
    ser = ser.rename('AVRG_BF')
    ser = ser.rolling(window=window, center=True).mean()
    ser = ser.fillna(0)

    index_list = range(0, len(bdamAtomList))
    index = pd.Series(index_list)
    index = index.rename('INDEX')

    df = pd.concat([df, ser, index], axis=1)

    df.loc[(df.AVRG_BF == 0) & (df.INDEX <= (math.floor(window/2)-1)),
           'AVRG_BF'] += df.BFAC.values[0:window].mean(axis=0)
    df.loc[(df.AVRG_BF == 0) & (df.INDEX >= (len(bdamAtomList) - math.floor(window/2))),
           'AVRG_BF'] += df.BFAC.values[(len(bdamAtomList)-window):len(bdamAtomList)].mean(axis=0)

    df = df.sort_values(by='ATMNUM', ascending=True)
    df = df.reset_index(drop=True)

    # The BDamage value of each atom in the DataFrame is calculated as the
    # ratio of its B factor value to its associated average B factor value.
    for index, atm in enumerate(bdamAtomList):
        atm.avrg_bf = df.AVRG_BF.values[index]
        atm.bd = atm.bFactor / atm.avrg_bf
