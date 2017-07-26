# RABDAM – identifying specific radiation damage in MX structures
A program to calculate the *B*<sub>Damage</sub> and *B*<sub>net</sub> metrics to quantify the extent of specific radiation damage present within an individual MX structure. Suitable for running on any standard format PDB file.

**NOTE:** These scripts are currently under development and are being updated regularly.


## Contents
-	[How to run in brief](#how-to-run-in-brief)
- [Background](#background)
-	[Usage](#usage)
    -	[System requirements](#system-requirements)
    -	[Running RABDAM from the command line](#running-rabdam-from-the-command-line)
    -	[Writing the RABDAM input file](#writing-the-rabdam-input-file)
- [An example RABDAM run](#an-example-rabdam-run)
-	[Queries](#queries)
-	[Contributors](#contributors)
-	[Citing RABDAM](#citing-rabdam)

## How to run in brief
RABDAM is a command line program. To run the program with its recommended default parameter values, execute:

`python rabdam.py –f XXXX`

, where XXXX is the 4 character PDB accession code of the MX structure under study. Alternatively, the user can provide RABDAM with a file path to a locally saved PDB file:

`python rabdam.py –f path/to/pdb_file.pdb`

See the *“Usage”* section below for further details.

## Background
During macromolecular crystallography (MX) data collection, X-rays are also absorbed by and deposit energy within the crystal under study, causing damage. This damage can result in localised chemical changes to the macromolecule copies, for instance to disulfide bond cleavage in proteins, *etc*. Such specific radiation damage manifestations can lead to incorrect biological conclusions being drawn from an MX structure if they are not identified and accounted for. Unfortunately, the high intensities of third generation synchrotron sources have resulted in specific radiation damage artefacts commonly being present in MX structures deposited in the Protein Data Bank (PDB) even at 100 K.

The chemical changes induced by specific radiation damage cause an accompanying increase in the atomic *B*<sub>factor</sub> values of affected sites. Multiple factors can affect an atom’s *B*<sub>factor</sub> value in addition to radiation damage however, the most important of which is its mobility: the increase in *B*<sub>factor</sub> caused by specific radiation damage is insufficiently large to distinguish damage from mobility.

There is a strong positive correlation between the mobility of an atom within a crystal structure and its packing density, i.e. the number of atoms present in its local environment. The *B*<sub>Damage</sub> metric is *B*<sub>factor</sub> corrected for packing density: specifically, the *B*<sub>Damage</sub> value for an atom *j* is calculated as the ratio between its *B*<sub>factor</sub> and the average *B*<sub>factor</sub> of atoms 1 to *n* which occupy a similar packing density environment to atom *j*.

The *B*<sub>Damage</sub> metric has been shown to identify expected sites of specific radiation damage in damaged datasets (Gerstel *et al.*, 2015).

RABDAM calculates the value of the *B*<sub>Damage</sub> metric for every atom within a standard format PDB file, according to the method summarised in the diagram below.


## Usage
#### System requirements
#### Running RABDAM from the command line
#### Writing the RABDAM input file

## An example RABDAM run

## Queries
Please email kathryn.shelley@queens.ox.ac.uk.

## Contributors
- Kathryn Shelley
- Tom Dixon
- Jonny Brooks-Bartlett

## Citing RABDAM
The initial development and testing of the BDamage metric is described in:

- Gerstel M, Deane CM, Garman EF (2015) Identifying and quantifying radiation damage at the atomic level. *J Synchrotron Radiat* **22**: 201-212. [http://dx.doi.org/doi:10.1107/S1600577515002131](http://dx.doi.org/doi:10.1107/S1600577515002131)

Please cite this paper if you use RABDAM to analyse specific radiation damage in your MX structure.

RABDAM is dependent upon the CCP4 suite program PDBCUR:

- Winn MD, Ballard CC, Cowtan KD, Dodson EJ, Emsley P, Evans PR, Keegan RM, Krissinel EB, Leslie AGW, McCoy A, McNicholas SJ, Murshudov GN, Pannu NS, Potterton EA, Powell HR, Read RJ, Vagin A, Wilson KS (2011) Overview of the CCP4 suite and current developments. *Acta Crystallogr Sect D Biol Crystallogr* **67**: 235-242.
