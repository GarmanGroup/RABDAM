# B_Damage
A program use to calculate B_Damage values for a biomolecular structure whose atomic coordinates are stored in PDB format.

** NOTE: The program is currently written in MATLAB but a version written in Python is currently under development.**

## Introduction
For molecular structures determined using X-ray crystallography each atom in the structure is assigned an ***atomic B factor*** value. This value effectively represents our level of uncertainty about the allocated position of that atom. The most mobile atoms will have the highest B factor values. However an atom's mobility could be due to different factors: increased thermal motion due to absorbed energy from the incident X-ray photons (a sign of *radiation damage*), or it could simply be due to the atom being in a highly flexible region of the protein and it is not surrounded by other 'stationary' atoms. ***B_Damage*** is a metric that attempts to deconvolute these factors to give an "effective" B factor value that is just a measure of the damage.   

**B_Damage** for a given atom is the ratio of the atom's B factor and the average B factor of atoms within a *similar packing density*. Atoms are grouped into **Similar packing density** environments, meaning that atoms in that group have a similar surrounding atomic environment which can be defined in several ways. The default here is that atoms with similar packing densities have a similar number of atoms within a given radius of the atom.

## Usage
**NOTE: This currently only refers to the version written in MATLAB. We'll update this when the Python version is complete.**    

***Need to complete this***

## Contributors
- Thomas Dixon
- Jonathan Brooks-Bartlett    

**Please cite** M. Gerstel, C. M. Deane and E.F. Garman. (2015). J. Synchrotron Radiation. **22**, 201-212 [http://dx.doi.org/doi:10.1107/S1600577515002131](http://dx.doi.org/doi:10.1107/S1600577515002131)
