# Copyright 2015 Thomas Dixon
# With thanks to Jonathan Brooks-Bartlett, Charles Bury, Markus Gerstel and Elspeth Garman

#Script to calculate B-damage for protein atoms
import sys
sys.path.insert(0,'./Subroutines')
from CalculateBdamage import cambda

cambda('2BN3')