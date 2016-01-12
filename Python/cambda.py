# Copyright 2015 Thomas Dixon
# With thanks to Jonathan Brooks-Bartlett, Charles Bury, Markus Gerstel and Elspeth Garman

#Script to calculate B-damage for protein atoms
from CalculateBdamage import cambda

#write your own inputs to this line, see the handbook for more details
cambda('4TOC', PDT=14, binSize=10, createAUCpdb=True, createTApdb=True)