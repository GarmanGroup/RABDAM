
# RABDAM
# Copyright (C) 2025 Garman Group, University of Oxford

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


def write_output_cif(atom_list, file_start, bdam):
    """
    Writes an mmCIF file for the input structure with an additional
    column of BDamage values for those atoms included in the calculation.
    """

    import copy

    cif_data = {
        'lineID': [None]*len(atom_list),
        'atomNum': [None]*len(atom_list),
        'origAtomType': [None]*len(atom_list),
        'conformer': [None]*len(atom_list),
        'origResiType': [None]*len(atom_list),
        'origChainID': [None]*len(atom_list),
        'entity_id': [None]*len(atom_list),
        'origResiNum': [None]*len(atom_list),
        'insCode': [None]*len(atom_list),
        'x': [None]*len(atom_list),
        'y': [None]*len(atom_list),
        'z': [None]*len(atom_list),
        'occupancy': [None]*len(atom_list),
        'bFactor': [None]*len(atom_list),
        'element': [None]*len(atom_list),
        'charge': [None]*len(atom_list),
        'resiNum': [None]*len(atom_list),
        'resiType': [None]*len(atom_list),
        'chainID': [None]*len(atom_list),
        'atomType': [None]*len(atom_list),
        'model_num': [None]*len(atom_list)
    }
    if bdam is True:
        cif_data['pd'] = [None]*len(atom_list)
        cif_data['avrg_bf'] = [None]*len(atom_list)
        cif_data['bd'] = [None]*len(atom_list)

    cif_column_widths = {
        'lineID': 0,
        'atomNum': 0,
        'origAtomType': 0,
        'conformer': 0,
        'origResiType': 0,
        'origChainID': 0,
        'entity_id': 0,
        'origResiNum': 0,
        'insCode': 0,
        'x': 0,
        'y': 0,
        'z': 0,
        'occupancy': 0,
        'bFactor': 0,
        'element': 0,
        'charge': 0,
        'resiNum': 0,
        'resiType': 0,
        'chainID': 0,
        'atomType': 0,
        'model_num': 0
    }
    if bdam is True:
        cif_column_widths['pd'] = 0
        cif_column_widths['avrg_bf'] = 0
        cif_column_widths['bd'] = 0

    for index, atom in enumerate(atom_list):
        cif_data['lineID'][index] = atom.lineID
        if len(cif_data['lineID'][index]) > cif_column_widths['lineID']:
            cif_column_widths['lineID'] = len(cif_data['lineID'][index])
        cif_data['atomNum'][index] = str(atom.atomNum)
        if len(cif_data['atomNum'][index]) > cif_column_widths['atomNum']:
            cif_column_widths['atomNum'] = len(cif_data['atomNum'][index])
        cif_data['origAtomType'][index] = atom.origAtomType
        if len(cif_data['origAtomType'][index]) > cif_column_widths['origAtomType']:
            cif_column_widths['origAtomType'] = len(cif_data['origAtomType'][index])
        cif_data['conformer'][index] = atom.conformer
        if len(cif_data['conformer'][index]) > cif_column_widths['conformer']:
            cif_column_widths['conformer'] = len(cif_data['conformer'][index])
        cif_data['origResiType'][index] = atom.origResiType
        if len(cif_data['origResiType'][index]) > cif_column_widths['origResiType']:
            cif_column_widths['origResiType'] = len(cif_data['origResiType'][index])
        cif_data['origChainID'][index] = atom.origChainID
        if len(cif_data['origChainID'][index]) > cif_column_widths['origChainID']:
            cif_column_widths['origChainID'] = len(cif_data['origChainID'][index])
        cif_data['entity_id'][index] = str(atom.entity_id)
        if len(cif_data['entity_id'][index]) > cif_column_widths['entity_id']:
            cif_column_widths['entity_id'] = len(cif_data['entity_id'][index])
        cif_data['origResiNum'][index] = atom.origResiNum
        if len(cif_data['origResiNum'][index]) > cif_column_widths['origResiNum']:
            cif_column_widths['origResiNum'] = len(cif_data['origResiNum'][index])
        cif_data['insCode'][index] = atom.insCode
        if len(cif_data['insCode'][index]) > cif_column_widths['insCode']:
            cif_column_widths['insCode'] = len(cif_data['insCode'][index])
        cif_data['x'][index] = str(atom.xyzCoords[0][0])
        if len(cif_data['x'][index]) > cif_column_widths['x']:
            cif_column_widths['x'] = len(cif_data['x'][index])
        cif_data['y'][index] = str(atom.xyzCoords[1][0])
        if len(cif_data['y'][index]) > cif_column_widths['y']:
            cif_column_widths['y'] = len(cif_data['y'][index])
        cif_data['z'][index] = str(atom.xyzCoords[2][0])
        if len(cif_data['z'][index]) > cif_column_widths['z']:
            cif_column_widths['z'] = len(cif_data['z'][index])
        cif_data['occupancy'][index] = str(atom.occupancy)
        if len(cif_data['occupancy'][index]) > cif_column_widths['occupancy']:
            cif_column_widths['occupancy'] = len(cif_data['occupancy'][index])
        cif_data['bFactor'][index] = str(atom.bFactor)
        if len(cif_data['bFactor'][index]) > cif_column_widths['bFactor']:
            cif_column_widths['bFactor'] = len(cif_data['bFactor'][index])
        cif_data['element'][index] = atom.element
        if len(cif_data['element'][index]) > cif_column_widths['element']:
            cif_column_widths['element'] = len(cif_data['element'][index])
        cif_data['charge'][index] = atom.charge
        if len(cif_data['charge'][index]) > cif_column_widths['charge']:
            cif_column_widths['charge'] = len(cif_data['charge'][index])
        cif_data['resiNum'][index] = str(atom.resiNum)
        if len(cif_data['resiNum'][index]) > cif_column_widths['resiNum']:
            cif_column_widths['resiNum'] = len(cif_data['resiNum'][index])
        cif_data['resiType'][index] = atom.resiType
        if len(cif_data['resiType'][index]) > cif_column_widths['resiType']:
            cif_column_widths['resiType'] = len(cif_data['resiType'][index])
        cif_data['chainID'][index] = atom.chainID
        if len(cif_data['chainID'][index]) > cif_column_widths['chainID']:
            cif_column_widths['chainID'] = len(cif_data['chainID'][index])
        cif_data['atomType'][index] = atom.atomType
        if len(cif_data['atomType'][index]) > cif_column_widths['atomType']:
            cif_column_widths['atomType'] = len(cif_data['atomType'][index])
        cif_data['model_num'][index] = str(atom.pdb_model_num)
        if len(cif_data['model_num'][index]) > cif_column_widths['model_num']:
            cif_column_widths['model_num'] = len(cif_data['model_num'][index])

        if bdam is True:
            cif_data['pd'][index] = str(int(atom.pd))
            if len(cif_data['pd'][index]) > cif_column_widths['pd']:
                cif_column_widths['pd'] = len(cif_data['pd'][index])
            cif_data['avrg_bf'][index] = str(round(atom.avrg_bf, 3))
            if len(cif_data['avrg_bf'][index]) > cif_column_widths['avrg_bf']:
                cif_column_widths['avrg_bf'] = len(cif_data['avrg_bf'][index])
            cif_data['bd'][index] = str(round(atom.bd, 3))
            if len(cif_data['bd'][index]) > cif_column_widths['bd']:
                cif_column_widths['bd'] = len(cif_data['bd'][index])

    for key, value_list in copy.deepcopy(cif_data).items():
        if set(value_list) == {''}:
            cif_data[key] = ['?']*len(atom_list)
        else:
            repl_value_list = []
            for val in value_list:
                if val == '':
                    repl_value_list.append('.'.ljust(cif_column_widths[key]))
                else:
                    repl_value_list.append(val.ljust(cif_column_widths[key]))
            cif_data[key] = repl_value_list

    # Writes output mmCIF file
    new_cif = open('%s.cif' % file_start, 'w')

    new_cif.write('data_{}\n'.format(file_start))
    new_cif.write('#\nloop_\n')
    new_cif.write(
        '_atom_site.group_PDB\n'
        '_atom_site.id\n'
        '_atom_site.type_symbol\n'
        '_atom_site.label_atom_id\n'
        '_atom_site.label_alt_id\n'
        '_atom_site.label_comp_id\n'
        '_atom_site.label_asym_id\n'
        '_atom_site.label_entity_id\n'
        '_atom_site.label_seq_id\n'
        '_atom_site.pdbx_PDB_ins_code\n'
        '_atom_site.Cartn_x\n'
        '_atom_site.Cartn_y\n'
        '_atom_site.Cartn_z\n'
        '_atom_site.occupancy\n'
        '_atom_site.B_iso_or_equiv\n'
    )
    if bdam is True:
        new_cif.write(
            '_atom_site.B_Damage\n'
            '_atom_site.packing_density\n'
            '_atom_site.average_B_factor\n'
        )
    new_cif.write(
        '_atom_site.pdbx_formal_charge\n'
        '_atom_site.auth_seq_id\n'
        '_atom_site.auth_comp_id\n'
        '_atom_site.auth_asym_id\n'
        '_atom_site.auth_atom_id\n'
        '_atom_site.pdbx_PDB_model_num\n'
    )

    for index in range(len(atom_list)):
        cif_line = [
            cif_data['lineID'][index],
            cif_data['atomNum'][index],
            cif_data['element'][index],
            cif_data['origAtomType'][index],
            cif_data['conformer'][index],
            cif_data['origResiType'][index],
            cif_data['origChainID'][index],
            cif_data['entity_id'][index],
            cif_data['origResiNum'][index],
            cif_data['insCode'][index],
            cif_data['x'][index],
            cif_data['y'][index],
            cif_data['z'][index],
            cif_data['occupancy'][index],
            cif_data['bFactor'][index]
        ]
        if bdam is True:
            cif_line += [
                cif_data['bd'][index],
                cif_data['pd'][index],
                cif_data['avrg_bf'][index]
            ]
        cif_line += [
            cif_data['charge'][index],
            cif_data['resiNum'][index],
            cif_data['resiType'][index],
            cif_data['chainID'][index],
            cif_data['atomType'][index],
            cif_data['model_num'][index]
        ]
        cif_line = ' '.join(cif_line)
        new_cif.write('{}\n'.format(cif_line))
    new_cif.write('#\n')

    new_cif.close()


def write_all_carbon_cif(atom_list, atom_id_list, file_start):
    """
    Writes an mmCIF file of carbon atoms at the xyz coordinates of each atom in
    'atom_list'. Suitable for unit cell, 3x3x3 unit cell and trimmed 3x3x3 unit
    cell only (to check that atoms are being translated, trimmed etc. correctly)
    """

    import copy

    cif_data = {
        'lineID': [None]*len(atom_list),
        'atomNum': [None]*len(atom_list),
        'origAtomType': [None]*len(atom_list),
        'conformer': [None]*len(atom_list),
        'origResiType': [None]*len(atom_list),
        'origChainID': [None]*len(atom_list),
        'entity_id': [None]*len(atom_list),
        'origResiNum': [None]*len(atom_list),
        'insCode': [None]*len(atom_list),
        'x': [None]*len(atom_list),
        'y': [None]*len(atom_list),
        'z': [None]*len(atom_list),
        'occupancy': [None]*len(atom_list),
        'bFactor': [None]*len(atom_list),
        'element': [None]*len(atom_list),
        'charge': [None]*len(atom_list),
        'resiNum': [None]*len(atom_list),
        'resiType': [None]*len(atom_list),
        'chainID': [None]*len(atom_list),
        'atomType': [None]*len(atom_list),
        'model_num': [None]*len(atom_list)
    }

    cif_column_widths = {
        'lineID': 0,
        'atomNum': 0,
        'origAtomType': 0,
        'conformer': 0,
        'origResiType': 0,
        'origChainID': 0,
        'entity_id': 0,
        'origResiNum': 0,
        'insCode': 0,
        'x': 0,
        'y': 0,
        'z': 0,
        'occupancy': 0,
        'bFactor': 0,
        'element': 0,
        'charge': 0,
        'resiNum': 0,
        'resiType': 0,
        'chainID': 0,
        'atomType': 0,
        'model_num': 0
    }

    for index, atom in enumerate(atom_list):
        cif_data['lineID'][index] = 'ATOM'
        cif_column_widths['lineID'] = len(cif_data['lineID'][index])
        cif_data['atomNum'][index] = str(atom_id_list[index])
        if len(cif_data['atomNum'][index]) > cif_column_widths['atomNum']:
            cif_column_widths['atomNum'] = len(cif_data['atomNum'][index])
        cif_data['origAtomType'][index] = 'C'
        cif_column_widths['origAtomType'] = len(cif_data['origAtomType'][index])
        cif_data['conformer'][index] = '?'
        cif_column_widths['conformer'] = len(cif_data['conformer'][index])
        cif_data['origResiType'][index] = 'GLY'
        cif_column_widths['origResiType'] = len(cif_data['origResiType'][index])
        cif_data['origChainID'][index] = 'A'
        cif_column_widths['origChainID'] = len(cif_data['origChainID'][index])
        cif_data['entity_id'][index] = '1'
        cif_column_widths['entity_id'] = len(cif_data['entity_id'][index])
        cif_data['origResiNum'][index] = '1'
        cif_column_widths['origResiNum'] = len(cif_data['origResiNum'][index])
        cif_data['insCode'][index] = '?'
        cif_column_widths['insCode'] = len(cif_data['insCode'][index])
        cif_data['x'][index] = str(atom_list[index][0])
        if len(cif_data['x'][index]) > cif_column_widths['x']:
            cif_column_widths['x'] = len(cif_data['x'][index])
        cif_data['y'][index] = str(atom_list[index][1])
        if len(cif_data['y'][index]) > cif_column_widths['y']:
            cif_column_widths['y'] = len(cif_data['y'][index])
        cif_data['z'][index] = str(atom_list[index][2])
        if len(cif_data['z'][index]) > cif_column_widths['z']:
            cif_column_widths['z'] = len(cif_data['z'][index])
        cif_data['occupancy'][index] = '1.00'
        cif_column_widths['occupancy'] = len(cif_data['occupancy'][index])
        cif_data['bFactor'][index] = '50.0'
        cif_column_widths['bFactor'] = len(cif_data['bFactor'][index])
        cif_data['element'][index] = 'C'
        cif_column_widths['element'] = len(cif_data['element'][index])
        cif_data['charge'][index] = '?'
        cif_column_widths['charge'] = len(cif_data['charge'][index])
        cif_data['resiNum'][index] = '1'
        cif_column_widths['resiNum'] = len(cif_data['resiNum'][index])
        cif_data['resiType'][index] = 'GLY'
        cif_column_widths['resiType'] = len(cif_data['resiType'][index])
        cif_data['chainID'][index] = 'A'
        cif_column_widths['chainID'] = len(cif_data['chainID'][index])
        cif_data['atomType'][index] = 'C'
        cif_column_widths['atomType'] = len(cif_data['atomType'][index])
        cif_data['model_num'][index] = '1'
        cif_column_widths['model_num'] = len(cif_data['model_num'][index])

    for key, value_list in copy.deepcopy(cif_data).items():
        if set(value_list) == {''}:
            cif_data[key] = ['?']*len(atom_list)
        else:
            repl_value_list = []
            for val in value_list:
                if val == '':
                    repl_value_list.append('.'.ljust(cif_column_widths[key]))
                else:
                    repl_value_list.append(val.ljust(cif_column_widths[key]))
            cif_data[key] = repl_value_list

    # Writes output mmCIF file
    new_cif = open('%s.cif' % file_start, 'w')

    new_cif.write('data_{}\n'.format(file_start))
    new_cif.write('#\nloop_\n')
    new_cif.write(
        '_atom_site.group_PDB\n'
        '_atom_site.id\n'
        '_atom_site.type_symbol\n'
        '_atom_site.label_atom_id\n'
        '_atom_site.label_alt_id\n'
        '_atom_site.label_comp_id\n'
        '_atom_site.label_asym_id\n'
        '_atom_site.label_entity_id\n'
        '_atom_site.label_seq_id\n'
        '_atom_site.pdbx_PDB_ins_code\n'
        '_atom_site.Cartn_x\n'
        '_atom_site.Cartn_y\n'
        '_atom_site.Cartn_z\n'
        '_atom_site.occupancy\n'
        '_atom_site.B_iso_or_equiv\n'
        '_atom_site.pdbx_formal_charge\n'
        '_atom_site.auth_seq_id\n'
        '_atom_site.auth_comp_id\n'
        '_atom_site.auth_asym_id\n'
        '_atom_site.auth_atom_id\n'
        '_atom_site.pdbx_PDB_model_num\n'
    )

    for index in range(len(atom_list)):
        cif_line = [
            cif_data['lineID'][index],
            cif_data['atomNum'][index],
            cif_data['element'][index],
            cif_data['origAtomType'][index],
            cif_data['conformer'][index],
            cif_data['origResiType'][index],
            cif_data['origChainID'][index],
            cif_data['entity_id'][index],
            cif_data['origResiNum'][index],
            cif_data['insCode'][index],
            cif_data['x'][index],
            cif_data['y'][index],
            cif_data['z'][index],
            cif_data['occupancy'][index],
            cif_data['bFactor'][index],
            cif_data['charge'][index],
            cif_data['resiNum'][index],
            cif_data['resiType'][index],
            cif_data['chainID'][index],
            cif_data['atomType'][index],
            cif_data['model_num'][index]
        ]
        cif_line = ' '.join(cif_line)
        new_cif.write('{}\n'.format(cif_line))
    new_cif.write('#\n')

    new_cif.close()


def calc_kde(bdam_vals):
    """
    Calculates a KDE using scipy
    """

    import numpy as np
    from scipy.stats import gaussian_kde

    bandwidth = gaussian_kde(dataset=bdam_vals).scotts_factor()
    tail_width = bandwidth*np.std(bdam_vals, ddof=0)
    x_min = min(bdam_vals) - (3*tail_width)
    x_max = max(bdam_vals) + (3*tail_width)

    kde = gaussian_kde(dataset=bdam_vals, bw_method='scott')
    x_values = np.linspace(x_min, x_max, 100)
    y_values = kde(x_values)

    return x_values, y_values


def calc_bnet(bnet_df, median):
    """
    Calculates the Bnet metric. To do this, a kernel density estimate (KDE) is
    calculated for the BDamage values of all glutamate and aspartate side chain
    carboxyl group oxygen atoms. The area under the KDE either side of the
    median BDamage values of all atoms in the structure is calculated; Bnet is
    then calculated as RHS / LHS, where RHS and LHS are the areas under the KDE
    to the right and left of the median, respectively.
    """
    
    x_values, y_values = calc_kde(bnet_df.BDAM.values)
    height = (x_values[-1]-x_values[0]) / (len(x_values)-1)

    total_area_LHS = 0
    for index, value in enumerate(y_values):
        if x_values[index] < median:
            area_LHS = (((y_values[index] + y_values[index+1]) / 2)
                        * height)
            total_area_LHS = total_area_LHS + area_LHS

    total_area_RHS = 0
    for index, value in enumerate(y_values):
        if x_values[index] >= median and index < len(y_values)-1:
            area_RHS = (((y_values[index] + y_values[index+1]) / 2)
                        * height)
            total_area_RHS = total_area_RHS + area_RHS

    # Calculates area ratio (= Bnet)
    bnet = total_area_RHS / total_area_LHS
    
    return bnet, x_values, y_values


def calc_bnet_percentile(bnet, resolution, bnet_database_df):
    """
    """

    import copy
    import numpy as np

    # Check that Bnet is not NaN or infinite, which will prevent Bnet_percentile
    # from being calculated
    if np.isnan(bnet) or np.isinf(bnet):
        print(
            'WARNING: Cannot calculate Bnet percentile - Bnet value is {}'.format(bnet)
        )
        bnet_percentile = np.nan
        return bnet_percentile
    
    # Check resolution has been correctly parsed from the input file
    if not type(resolution) in [int, float, np.float32, np.float64]:
        print(
            'WARNING: Cannot calculate Bnet percentile - resolution value is {}'
            ', type {}, RABDAM expects a numeric value'.format(resolution, type(resolution))
        )
        bnet_percentile = np.nan
        return bnet_percentile
    

    # Find the 1000 structures closest in resolution to the input structure
    res_array = copy.deepcopy(bnet_database_df['Resolution (A)']).to_numpy()
    surr_struct_indices = []
    surr_struct_vals = []
    for num in range(1000):
        nearest_res_index = (np.abs(res_array-resolution)).argmin()
        nearest_resolution = res_array[nearest_res_index]
        surr_struct_indices.append(nearest_res_index)
        surr_struct_vals.append(nearest_resolution)
        res_array[nearest_res_index] = np.inf

    # Include any other structures that fall in this resolution range
    min_val = min(surr_struct_vals)
    max_val = max(surr_struct_vals)
    for index, num in np.ndenumerate(res_array):
        if num == min_val or num == max_val:
            surr_struct_indices.append(index[0])

    # Calculate Bnet_percentile
    surr_struct_df = copy.deepcopy(bnet_database_df).iloc[surr_struct_indices].reset_index(drop=True)
    bnet_array = np.sort(surr_struct_df['Bnet'].to_numpy())
    num_entries = bnet_array.shape[0]
    nearest_bnet_index = (np.abs(bnet_array-bnet)).argmin()
    bnet_percentile = ((nearest_bnet_index + 1) / num_entries) * 100

    return bnet_percentile


class generate_output_files(object):
    def __init__(self, out_file_start, pdb_code, df, prot_or_na, resolution):
        self.out_file_start = out_file_start
        self.pdb_code = pdb_code
        self.df = df
        self.prot_or_na = prot_or_na
        self.resolution = resolution

    def make_csv(self, window):
        """
        Returns a csv file containing a complete set of atom information
        (including both that provided in the input PDB / mmCif file and also
        the BDamage values calculated by RABDAM) for all atoms considered for
        BDamage analysis. (This provides the user with a copy of the raw data
        which they can manipulate as they wish.)
        """

        newFile = open('%s_BDamage.csv' % self.out_file_start, 'w')

        # Defines column header abbreviations
        newFile.write('REC = RECORD NAME\n'
                      'ATMNUM = ATOM SERIAL NUMBER\n'
                      'ATMNAME = ATOM NAME\n'
                      'CONFORMER = ALTERNATE LOCATION INDICATOR\n'
                      'RESNAME = RESIDUE NAME\n'
                      'CHAIN = CHAIN IDENTIFIER\n'
                      'RESNUM = RESIDUE SEQUENCE NUMBER\n'
                      'INSCODE = INSERTION CODE'
                      'XPOS = ORTHOGONAL COORDINATES FOR X IN ANGSTROMS\n'
                      'YPOS = ORTHOGONAL COORDINATES FOR Y IN ANGSTROMS\n'
                      'ZPOS = ORTHOGONAL COORDINATES FOR Z IN ANGSTROMS\n'
                      'OCC = OCCUPANCY\n'
                      'BFAC = B FACTOR (TEMPERATURE FACTOR)\n'
                      'ELEMENT = ELEMENT SYMBOL\n'
                      'CHARGE = CHARGE ON ATOM\n'
                      'PD = PACKING DENSITY (ATOMIC CONTACT NUMBER)\n')
        newFile.write('AVRG_BF = AVERAGE B FACTOR FOR ATOMS IN A SIMILAR '
                      'PACKING DENSITY ENVIRONMENT (SLIDING WINDOW SIZE '
                      '= %s)\n' % window)
        newFile.write('BDAM = B DAMAGE VALUE\n'
                      '\n')

        # Writes properties of each atom considered for BDamage analysis.
        newFile.write(self.df.to_csv(index=False))
        newFile.close()

    def gen_kde_plot(self, highlightAtoms):
        """
        Returns a kernel density estimate of the BDamage values of every atom
        considered for BDamage analysis. Any atom whose number is listed
        in the highlightAtoms option in the input file will be marked on the
        plot. (Note that it is recommended no more than 6 atoms are listed
        in the highlightAtoms option in the input file (beyond 6 atoms, the
        colour scheme will repeat itself, and in addition the key may not fit
        onto the graph).)
        """

        import matplotlib.pyplot as plt

        plt.clf()
        # Calculate KDE
        x_values, y_values = calc_kde(self.df.BDAM.values)
        # Plot rugplot
        plt.plot(
            self.df.BDAM.values, [0.05]*len(self.df.BDAM.values), '|',
            color='#212f3d', alpha=0.25
        )
        # Plot KDE
        plt.plot(x_values, y_values, color='#212f3d')

        # Marks on the positions of any atoms whose numbers are listed in the
        # highlightAtoms option specified in the input file.
        highlighted_atoms = []
        for number in highlightAtoms:
            index = self.df.ATMNUM.tolist().index(int(number))
            b_dam_value = self.df.BDAM.values[index]
            line, = plt.plot([b_dam_value, b_dam_value], [0, max(y_values)],
                             linewidth=2, label=' atom ' + str(number) +
                             '\n BDamage = {:.2f}'.format(b_dam_value))
            highlighted_atoms.append(line)

        if len(highlighted_atoms) >= 1:
            plt.legend(handles=highlighted_atoms)
        plt.xlabel('B Damage')
        plt.ylabel('Normalised Frequency')
        plt.title(self.pdb_code + ' BDamage kernel density plot')
        plt.grid(visible=True, color='white')
        ax = plt.gca()
        ax.set_facecolor('#f0f3f7')
        plt.savefig(self.out_file_start + '_BDamage.svg')

    def gen_Bnet_out_files(self, bnet_df, median, prot_or_na, bnet_database_df):
        """
        """

        import os
        import matplotlib.pyplot as plt
        import pandas as pd

        # Calculate Bnet
        bnet, x_values, y_values = calc_bnet(bnet_df, median)
        # Calculate Bnet_percentile
        bnet_percentile = calc_bnet_percentile(
            bnet, self.resolution, bnet_database_df
        )

        # Plot KDE and rugplot
        plt.clf()
        plt.plot(
            bnet_df.BDAM.values, [0.05]*len(bnet_df.BDAM.values), '|',
            color='#212f3d', alpha=0.5
        )
        plt.plot(x_values, y_values, color='#212f3d')
        plt.xlabel('B Damage')
        plt.ylabel('Normalised Frequency')
        plt.title(self.pdb_code + ' Bnet kernel density plot')

        # Annotate KDE plot with Bnet and Bnet_percentile values
        plt.annotate('Bnet = {:.1f}'.format(bnet),
                     xy=(max(x_values)*0.65, max(y_values)*0.9),
                     fontsize=10)
        plt.annotate('Bnet_percentile = {:.1f}'.format(bnet_percentile),
                     xy=(max(x_values)*0.65, max(y_values)*0.85),
                     fontsize=10)
        plt.annotate('Median = {:.2f}'.format(median),
                     xy=(max(x_values)*0.65, max(y_values)*0.8),
                     fontsize=10)
        plt.grid(visible=True, color='white')
        ax = plt.gca()
        ax.set_facecolor('#f0f3f7')
        plt.savefig(self.out_file_start + '_Bnet_{}.svg'.format(prot_or_na))

        # Save output scores as a csv file and as a pickled dataframe
        if not os.path.isfile('Logfiles/Bnet_{}.csv'.format(prot_or_na)):
            bnet_list = open('Logfiles/Bnet_{}.csv'.format(prot_or_na), 'w')
            bnet_list.write('PDB' + ',')
            bnet_list.write('Bnet' + ',')
            bnet_list.write('Bnet percentile' + '\n')
            bnet_list.close()
        bnet_list = open('Logfiles/Bnet_{}.csv'.format(prot_or_na), 'a')
        bnet_list.write('%s,%s,%s\n' % (self.pdb_code, bnet, bnet_percentile))
        bnet_list.close()

        out_df = pd.DataFrame({'PDB': [self.pdb_code],
                               'Bnet': [bnet],
                               'Bnet_percentile': [bnet_percentile]})
        if os.path.isfile('Logfiles/Bnet_{}.pkl'.format(prot_or_na)):
            old_out_df = pd.read_pickle('Logfiles/Bnet_{}.pkl'.format(prot_or_na))
        else:
            old_out_df = pd.DataFrame({'PDB': [],
                                       'Bnet': [],
                                       'Bnet_percentile': []})
        new_out_df = pd.concat(
            [old_out_df, out_df], axis=0, ignore_index=True
        ).reset_index(drop=True)
        new_out_df.to_pickle('Logfiles/Bnet_{}.pkl'.format(prot_or_na))

        return bnet, bnet_percentile

    def calculate_Bnet_and_Bnet_percentile(self, bnet_database_df):
        """
        Plots a kernel density estimate of the BDamage values of Glu O and
        Asp O atoms. The summary metric Bnet is then calculated as the ratio
        of the areas under the curve either side of the median (of the
        overall BDamage distribution).
        """

        import sys
        # Ensures that seaborn uses scipy for kde bandwidth calculation rather
        # than statsmodels => must be run before import seaborn
        sys.modules['statsmodels']=None
        import pandas as pd
        import matplotlib.pyplot as plt

        # Selects Glu / Asp (both the L- and D-amino acids) terminal oxygen
        # atoms from complete DataFrame.
        if self.prot_or_na in ['protein', 'proteinna']:
            prot_df = self.df[(self.df.RESNAME.isin(['GLU', 'ASP', 'DGL', 'DAS']))
                            & (self.df.ATMNAME.isin(['OE1', 'OE2', 'OD1', 'OD2']))]
        else:
            prot_df = pd.DataFrame({})
        if self.prot_or_na in ['nucleicacid', 'na', 'proteinna']:
            # Selects atoms of sugar-phosphate C-O bonds from complete DataFrame.
            na_df = self.df[self.df.ATMNAME.isin(["O3'", "O5'", "C3'", "C5'"])]
        else:
            na_df = pd.DataFrame({})

        # Calculates median of BDamage distribution
        median = self.df.BDAM.median()

        if prot_df.empty and na_df.empty:
            print('\n\nERROR: No sites used for Bnet calculation present in structure\n')
            return

        prot_bnet = None
        prot_bnet_percentile = None
        if not prot_df.empty:
            prot_bnet, prot_bnet_percentile = self.gen_Bnet_out_files(
                prot_df, median, 'protein', bnet_database_df
            )

        na_bnet = None
        na_bnet_percentile = None
        if not na_df.empty:
            na_bnet, na_bnet_percentile = self.gen_Bnet_out_files(
                na_df, median, 'NA', bnet_database_df
            )

        return prot_bnet, prot_bnet_percentile, na_bnet, na_bnet_percentile

    def write_html_summary(self, highlightAtoms):
        """
        Writes an html file summarising the results generated by RABDAM.
        """

        import os
        import requests
        import time
        import pandas as pd
        import numpy as np

        # Records time (for reporting in summary html file)
        tm = time.localtime()
        day = tm.tm_mday
        month = tm.tm_mon
        year = tm.tm_year
        hour = tm.tm_hour
        mins = tm.tm_min
        secs = tm.tm_sec

        # Filters the complete DataFrame to retain only the atoms specified to
        # be highlighted (default = none) by the user in the input file
        if len(highlightAtoms) != 0:
            highlightAtoms = [int(number) for number in highlightAtoms]
            sub_df_highlight = self.df[self.df.ATMNUM.isin(highlightAtoms)]
            sub_df_highlight = sub_df_highlight.round({'AVRG BF': 2, 'BDAM': 2})

        # Filter the complete DataFrame to retain the atoms with highest
        # BDamage values (> 1.96o = 2-tailed 95% confidence interval)
        sorted_df = self.df.sort_values(by='BDAM', ascending=False)
        b_dam_values = sorted_df['BDAM']
        ln_b_dam_values = np.log(b_dam_values)
        mean = np.mean(ln_b_dam_values)
        std_dev = np.std(ln_b_dam_values)
        cut = mean + (1.96*std_dev)
        cut_index = 0
        for index, value in enumerate(ln_b_dam_values.tolist()):
            if value > cut:
                cut_index = index
            else:
                break
        sub_df_top_site = sorted_df.drop(sorted_df.index[cut_index+1:])
        sub_df_top_site = sub_df_top_site.round({'AVRG BF': 2, 'BDAM': 2})

        # Filters the complete DataFrame to retain only Glu and Asp (L- and
        # D-isomers) terminal oxygen atoms.
        sub_df_prot = self.df[(self.df.RESNAME.isin(['GLU', 'ASP', 'DGL', 'DAS']))
                              & (self.df.ATMNAME.isin(['OE1', 'OE2', 'OD1', 'OD2']))]
        sub_df_prot = sub_df_prot.round({'AVRG BF': 2, 'BDAM': 2})
        # Filters the complete DataFrame to retain only atoms of
        # sugar-phosphate C-O bonds
        sub_df_na = self.df[self.df.ATMNAME.isin(["O3'", "O5'", "C3'", "C5'"])]
        sub_df_na = sub_df_na.round({'AVRG BF': 2, 'BDAM': 2})

        # Opens summary html file
        html_file = open(self.out_file_start + '_BDamage.html', 'w')
    
        # Try first to load html files from the internet
        js_file = 'https://raw.githubusercontent.com/GarmanGroup/RABDAM/master/rabdam/Subroutines/HTML_stylesheet.js'
        js_header = requests.get(js_file)
        if js_header.status_code == 200:
            js_text = js_header.text
        else:
            js_file = 'file://%s/HTML_stylesheet.js' % os.path.dirname(os.path.abspath(__file__))
            with open(js_file, 'r') as f:
                js_text = f.read()

        css_file = 'https://raw.githubusercontent.com/GarmanGroup/RABDAM/master/rabdam/Subroutines/HTML_stylesheet.css'
        css_header = requests.get(css_file)
        if css_header.status_code == 200:
            css_text = css_header.text
        else:
            css_file = 'file://%s/HTML_stylesheet.css' % os.path.dirname(os.path.abspath(__file__))
            with open(css_file, 'r') as f:
                css_text = f.read()

        jquery_file = 'http://ajax.googleapis.com/ajax/libs/jquery/1.10.2/jquery.min.js'
        jquery_header = requests.get(jquery_file)
        if jquery_header.status_code == 200:
            jquery_text = jquery_header.text
        else:
            jquery_file = 'file://%s/jQuery_stylesheet.js' % os.path.dirname(os.path.abspath(__file__))
            with open(jquery_file, 'r') as f:
                jquery_text = f.read()

        # Specifies html file header information
        html_file.write('<!DOCTYPE html>\n'
                        '<html>\n'
                        '  <head>\n'
                        '    <script type="text/javascript"> %s </script>\n' % jquery_text)
        html_file.write('    <script type="text/javascript"> %s </script>\n' % js_text)
        html_file.write('    <style> %s </style>\n' % css_text)
        html_file.write('    <title>'+self.pdb_code.replace('_', ' ')+' BDamage summary file</title>\n'
                        '  </head>\n')
        # Writes html file title
        html_file.write('  <body>\n'
                        '    <div class="file_name">\n'
                        '      <h1>'+self.pdb_code.replace('_', ' ')+' B<sub>Damage</sub> summary file</h1>\n')
        html_file.write('      <p id="file_info">Created on %02.0f/%02.0f/%04.0f at %02.0f:%02.0f:%02.0f.</p>\n'
                        % (day, month, year, hour, mins, secs))
        html_file.write('      <p id="file_info">**NOTE: Table header abbreviations are as defined in '
                        '<a href="%s_BDamage.csv">%s_Bdamage.csv</a>.**</p>\n' % (self.pdb_code, self.pdb_code))
        html_file.write('    </div>\n')
        # Writes html file BDamage summary
        html_file.write('    <div>\n'
                        '      <ul class="main_list">\n'
                        '        <li><h2>B<sub>Damage</sub> distribution</h2></li>\n'
                        '        <li><img class="graphs" src="%s_BDamage.svg" /></li>\n' % self.pdb_code)
        if len(highlightAtoms) != 0:
            html_file.write('        <div>\n'
                            '          <ul class="sublist">\n'
                            '            <li><h3>Properties of highlighted atoms</h3></li>\n'
                            '            <li>'+sub_df_highlight.to_html(index=False)+'</li></br>\n'
                            '          </ul>\n'
                            '        </div>\n')
        html_file.write('        <div>\n'
                        '          <ul class="sublist">\n'
                        '            <li><h3>Properties of sites with B<sub>Damage</sub> values larger than the 95% confidence limit (2-tailed) threshold</h3></li>\n'
                        '            <li>'+sub_df_top_site.to_html(index=False)+'</li>\n'
                        '          </ul>\n'
                        '        </div>\n'
                        '      </ul>\n'
                        '    </div>\n')
        # Writes html file Bnet (protein) summary
        if os.path.isfile('%s_Bnet_Protein.svg' % self.out_file_start):
            html_file.write('    <div>\n'
                            '      <ul class="main_list">\n'
                            '        <li><h2>B<sub>net</sub> distribution (protein)</h2></li>\n'
                            '        <li><img class="graphs" src="%s_Bnet_Protein.svg" /></li>\n' % self.pdb_code)
            html_file.write('        <ul class="sublist">\n'
                            '          <li><h3>Properties of atoms used in calculation of B<sub>net</sub> (protein) metric</h3></li>\n'
                            '          <li>'+sub_df_prot.to_html(index=False)+'</li>\n'
                            '        </ul>\n'
                            '      </ul>\n'
                            '    </div>\n')
        # Writes html file Bnet (nucleic acid) summary
        if os.path.isfile('%s_Bnet_NA.svg' % self.out_file_start):
            html_file.write('    <div>\n'
                            '      <ul class="main_list">\n'
                            '        <li><h2>B<sub>net</sub> distribution (nucleic acid)</h2></li>\n'
                            '        <li><img class="graphs" src="%s_Bnet_NA.svg" /></li>\n' % self.pdb_code)
            html_file.write('        <ul class="sublist">\n'
                            '          <li><h3>Properties of atoms used in calculation of B<sub>net</sub> (nucleic acid) metric</h3></li>\n'
                            '          <li>'+sub_df_na.to_html(index=False)+'</li>\n'
                            '        </ul>\n'
                            '      </ul>\n'
                            '    </div>\n')
        # Concludes and closes summary html file
        html_file.write('  </body>\n'
                        '</html>\n')
        html_file.close()
