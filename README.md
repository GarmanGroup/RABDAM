# RABDAM – identifying specific radiation damage in MX structures
A program to calculate the *B*<sub>Damage</sub> and *B*<sub>net</sub> metrics to quantify the extent of specific radiation damage present within an individual MX structure. Suitable for running on any standard format PDB file.

\*\***NOTE:** These scripts are under development, and are updated regularly. The program is currently being extended to incorporate nucleic acids analysis. Whilst these new capabilities are being tested, it is strongly recommended that presently RABDAM is only used to assess damage to protein crystal structures.\*\*


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

, where XXXX is the 4 character PDB accession code of the MX structure under study. Alternatively, you can provide RABDAM with a file path to a locally saved PDB file:

`python rabdam.py –f path/to/pdb_file.pdb`

See the “*Usage*” section below for further details.

## Background
During macromolecular crystallography (MX) data collection, X-rays are also absorbed by and deposit energy within the crystal under study, causing damage. This damage can result in localised chemical changes to the crystalline macromolecule copies, for example to disulfide bond cleavage in proteins, *etc*. Such *specific radiation damage* manifestations can lead to incorrect biological conclusions being drawn from an MX structure if they are not identified and accounted for. Unfortunately, the high intensities of third generation synchrotron sources have resulted in specific radiation damage artefacts commonly being present in MX structures deposited in the Protein Data Bank (PDB) even at 100 K.

The chemical changes induced by specific radiation damage cause an accompanying increase in the atomic *B*<sub>factor</sub> values of affected sites. Multiple factors can affect an atom’s *B*<sub>factor</sub> value in addition to radiation damage however, the most important of these being its mobility. The increase in *B*<sub>factor</sub> caused by specific radiation damage is insufficiently large to distinguish damage from mobility.

There is a strong positive correlation between the mobility of an atom within a crystal structure and its packing density, *i.e.* the number of atoms present in its local environment. The *B*<sub>Damage</sub> metric is *B*<sub>factor</sub> corrected for packing density: specifically, the *B*<sub>Damage</sub> value of an atom *j* is equal to the ratio of its *B*<sub>factor</sub> to the average *B*<sub>factor</sub> of atoms 1 to *n* which occupy a similar packing density environment to atom *j*. The *B*<sub>Damage</sub> metric has been shown to identify expected sites of specific radiation damage in damaged MX structures (Gerstel *et al.*, 2015).

![Images/BDamage_equation.png](Images/BDamage_equation.png)

The method of calculating an atom’s *B*<sub>Damage</sub> value is summarised in the diagram below:

___

![Images/BDamage_methodology.png](Images/BDamage_methodology.png)
Calculation of the *B*<sub>Damage</sub> metric. From an input PDB file of the asymmetric unit of a macromolecule of interest, RABDAM **(A)** generates a copy of the unit cell, followed by **(B)** a 3x3x3 assembly of unit cells. **(C)** Atoms in the 3x3x3 unit cell assembly that lie further than 14 Å from the asymmetric unit are discounted. **(D)** The packing density of an atom *j* in the asymmetric unit is calculated as the number of atoms within a 14 Å radius. **(E)** Asymmetric unit atoms are ordered by packing density; the *B*<sub>Damage</sub> value of atom *j* is then calculated as the ratio of its *B*<sub>factor</sub> to the average of the *B*<sub>factor</sub> values of atoms grouped, via sliding window, as occupying a similar packing density environment. (Hydrogen atoms are not considered in the calculation of *B*<sub>Damage</sub>.)

___

The *B*<sub>net</sub> metric is a derivative of the (per-atom) *B*<sub>Damage</sub> metric that summarises in a single value the overall extent of specific radiation damage suffered by an MX structure. One of the best-characterised chemical changes resulting from specific radiation damage that occurs within proteins\* is the decarboxylation of Glu and Asp residues; the *B*<sub>net</sub> metric is calculated from a kernel density plot of the *B*<sub>Damage</sub> values of a structure’s Glu and Asp terminal oxygen atoms as the ratio of the area either side of the median of the (overall) *B*<sub>Damage</sub> distribution.

(\* An equivalent of this protein-specific *B*<sub>net</sub> metric for nucleic acids is currently being developed - see program description.)

The method of calculating the *B*<sub>net</sub> value for a protein structure is summarised in the diagram below:

___

![Images/Bnet_calculation.png](Images/Bnet_calculation.png)

___

RABDAM calculates the values of the *B*<sub>Damage</sub> and *B*<sub>net</sub> metrics for a standard format PDB file, as detailed in the following sections.

## Usage
#### System requirements
___
#### Running RABDAM from the command line
___
#### Writing the RABDAM input file
___

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
