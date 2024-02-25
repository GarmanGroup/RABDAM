
# RABDAM
# Copyright (C) 2024 Garman Group, University of Oxford

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


from phenix.program_template import ProgramTemplate
import libtbx


# Define program citations
program_citations = libtbx.phil.parse('''
    citation {
    article_id = RABDAM
    authors = Shelley KL, Dixon TPE, Brooks-Bartlett JC, Garman EF
    title = RABDAM: quantifying specific radiation damage in individual protein crystal structures.
    journal = Journal of Applied Crystallography
    volume = 51
    pages = 552-559
    year = 2018
    doi_id = "10.1107/S1600576718002509"
    pmid = 29657569
    external = True
    }
    citation {
    article_id = Bnet and Bnet-percentile
    authors = Shelley KL, Garman EF
    title = Quantifying and comparing radiation damage in the Protein Data Bank.
    journal = Nature Communications
    volume = 13
    pages = 1314
    year = 2022
    doi_id = "10.1038/s41467-022-28934-0"
    pmid = 35288575
    external = True
    }
    citation {
    article_id = BDamage
    authors = Gerstel M, Deane CM, Garman EF
    title = Identifying and quantifying radiation damage at the atomic level.
    journal = Journal of Synchrotron Radiation
    volume = 22
    pages = 201-212
    year = 2015
    doi_id = "10.1107/S1600577515002131"
    pmid =  25723922
    external = True
    }
''')


# Define the input parameters
# Only the first two parameters are relevant for a GUI - default values will be
# acceptable for vast majority of cases
master_phil_str = '''
    input_path = None
        .type = path
    output_dir = None
        .type = path
    batch_run = True
        .type = bool
    overwrite = True
        .type = bool
    output_files = "all"
        .type = str
    filter = True
        .type = bool
    temperature = None
        .type = str
    resolution = None
        .type = float
    pdt = 7.0
        .type = float
    window = 0.02
        .type = float
    hetatm = False
        .type = bool
    remove_atoms_list = "[]"
        .type = str
    add_atoms_list = "[]"
        .type = str
    highlight_atoms_list = "[]"
        .type = str
    save_orig_pdb = False
        .type = bool
    save_au_pdb = False
        .type = bool
    save_uc_pdb = False
        .type = bool
    save_auc_pdb = False
        .type = bool
    save_ta_pdb = False
        .type = bool
    test = False
        .type = bool
    run_type = "full"
        .type = str
    phenix_import = True
        .type = bool
'''


def convert_input_str_to_list(input_string, variable_name):
    """
    Converts list in string format ("[1, 2, 3]") into list ([1, 2, 3])
    """

    input_string = str(input_string).strip()
    output_list = []

    if input_string not in ['', '[]']:
        input_string = input_string.lstrip('[').rstrip(']')
        input_list = input_string.split(',')  # Necessary to use "," rather 
        # than ";" as phenix throws an error if include ";" in the string
        for item in input_list:
            num_range = item.split('-').strip()
            if len(num_range) == 2:
                try:
                    min_val = int(num_range[0])
                    max_val = int(num_range[1])
                    full_num_range = range(min_val, (max_val+1))
                    for number in full_num_range:
                        output_list.append(str(number))
                except ValueError:
                    raise ValueError(
                        'Unrecognised input: {} for {} - if input contains "-",'
                        ' expecting a numeric range'.format(num_range, variable_name)
                    )
            elif len(num_range) == 1:
                output_list.append(str(num_range[0]))
            else:
                raise ValueError(
                    'Unrecognised input: {} for {} - if input contains "-", '
                    'expecting a numeric range'.format(num_range, variable_name)
                )

    return output_list


def convert_input_temp_to_float(temp_str):
    """
    Converts input temperature into a float value
    """

    temp_float = None
    temp_str = str(temp_str).lower().strip()

    if temp_str == 'cryo':
        temp_float = 100
    elif temp_str == 'none':
        temp_float = None
    else:
        try:
            temp_float = float(temp_str.rstrip('K'))
        except ValueError:
            raise ValueError(
                'Value provided for temperature unrecognised: {}\n'
                'Expect to be set to "cryo", or to a temperature value in '
                'Kelvin'.format(temp_str)
            )

    return temp_float


class Program(ProgramTemplate):

    # Program description
    description = (
        "RABDAM: A program to help users to assess protein crystal structures "
        "for radiation damage by calculating the BDamage, Bnet and "
        "Bnet-percentile metrics.\nBecause these metrics are Bfactor-based, "
        "they should only be calculated for models that have been sufficiently "
        "well-refined that they either have already been deposited in the PDB, "
        "or are close to being ready for PDB deposition.\nThe Bnet and "
        "Bnet-percentile metrics have currently been validated for "
        "cryo-temperature (80 K - 120 K) protein crystal structures only.\n"
    )

    # Define the data types that will be used
    datatypes = ["phil"]

    # Define program input parameters
    master_phil_str = master_phil_str

    # Define program citations
    citations = program_citations

    # Define how to determine if inputs are ok
    def validate(self):
        # Expect input_path to exist or to be a PDB accession code, and
        # output_dir to exist

        import os
        from libtbx.utils import Sorry

        input_path_error = False
        output_dir_error = False

        if self.params.input_path is None:
            input_path_error = True
        else:
            if len(self.params.input_path) != 4:  # Allows RABDAM to be run with
            # a PDB accession code
                if not os.path.isfile(self.params.input_path):
                    input_path_error = True

        if self.params.output_dir is None:
            output_dir_error = True
        else:
            if not os.path.isdir(self.params.output_dir):
                output_dir_error = True

        if input_path_error is True:
            raise Sorry(
                'Path to input model is not recognised:\n'
                '{}'.format(self.params.input_path)
            )
        if output_dir_error is True:
            raise Sorry(
                'Path to output directory is not recognised:\n'
                '{}'.format(self.params.output_dir)
            )

    # Run the program
    def run(self):
        from phenix.rabdam.Subroutines.CalculateBDamage import run_rabdam

        # Convert command line strings into lists where appropriate
        self.params.remove_atoms_list = convert_input_str_to_list(
            self.params.remove_atoms_list, "remove_atoms_list"
        )
        self.params.add_atoms_list = convert_input_str_to_list(
            self.params.add_atoms_list, "add_atoms_list"
        )
        self.params.highlight_atoms_list = convert_input_str_to_list(
            self.params.highlight_atoms_list, "highlight_atoms_list"
        )
        # Convert command line temperature from string to float/None
        self.params.temperature = convert_input_temp_to_float(
            self.params.temperature
        )

        # Initialises rabdam object
        rabdam_obj = run_rabdam(
            pathToInput=self.params.input_path,
            outputDir=self.params.output_dir,
            batchRun=self.params.batch_run,
            overwrite=self.params.overwrite,
            outFiles=self.params.output_files,
            filterInput=self.params.filter,
            temperature=self.params.temperature,
            resolution=self.params.resolution,
            PDT=self.params.pdt,
            windowSize=self.params.window,
            HETATM=self.params.hetatm,
            removeAtoms=self.params.remove_atoms_list,
            addAtoms=self.params.add_atoms_list,
            highlightAtoms=self.params.highlight_atoms_list,
            createOrigpdb=self.params.save_orig_pdb,
            createAUpdb=self.params.save_au_pdb,
            createUCpdb=self.params.save_uc_pdb,
            createAUCpdb=self.params.save_auc_pdb,
            createTApdb=self.params.save_ta_pdb,
            phenixImport=self.params.phenix_import
        )

        # Runs RABDAM to calculate BDamage and Bnet values
        success = True
        bnet = None
        bnet_percentile = None
        if self.params.run_type in ["full", "dataframe"]:
            success = rabdam_obj.rabdam_dataframe(test=self.params.test)

        if success is False:
            print("RABDAM failed to run for {}".format(self.params.input_path))
        else:
            if self.params.run_type in ["full", "analysis"]:
                # For now RABDAM only considers protein Bnet and Bnet percentile
                # values, so the nucleic acid Bnet and Bnet percentile values
                # returned are ignored
                bnet, bnet_percentile, _, _ = rabdam_obj.rabdam_analysis()
        self.result = {
            "bnet": bnet,
            "bnet_percentile": bnet_percentile,
        }

    # Return Bnet and Bnet_percentile values
    def get_results(self):
        from libtbx import group_args

        results_obj = group_args(
            bnet=self.result["bnet"],
            bnet_percentile=self.result["bnet_percentile"]
        )
        return results_obj


# Runs RABDAM from the "Program" class
def main():
    from iotbx.cli_parser import run_program
    run_program(program_class=Program)


# Runs "main" function if rabdam_phenix_template.py is run as a script 
if __name__ == "__main__":
    main()