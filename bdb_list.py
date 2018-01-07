
# RABDAM
# Copyright (C) 2017 Garman Group, University of Oxford

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

# A list of PDB accession codes, downloaded from the B-factor Databank
# (http://www.cmbi.umcn.nl/bdb/) and the RSCB PDB website
# (https://www.rcsb.org/pdb/), of structures refined with full isotropic
# atomic B-factors, and hence are suitable for RABDAM analysis.

def rabdam_compatible_structures():
    compat_struct_list = []
    return compat_struct_list
