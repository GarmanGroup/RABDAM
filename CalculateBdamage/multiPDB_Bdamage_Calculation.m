%% Script to calculate Bdamage for several pdb entries

clear all
close all
clc

%% Choose pdb entries

pdbEntries.pdb1 = '2bn3';
pdbEntries.pdb2 = '2bn1';
pdbEntries.pdb3 = '1f34';
pdbEntries.pdb4 = '2jhu';
pdbEntries.pdb5 = '4cyn';
pdbEntries.pdb6 = '2p1h';
pdbEntries.pdb7 = '4to2';
pdbEntries.pdb8 = '2blo';
pdbEntries.pdb9 = '2blq';
pdbEntries.pdb10 = '1gd6';
pdbEntries.pdb11 = '2blx';
pdbEntries.pdb12 = '2bly';
pdbEntries.pdb13 = '2blp';
pdbEntries.pdb14 = '2blz';
pdbEntries.pdb15 = '2blr';
pdbEntries.pdb16 = '2blu';
pdbEntries.pdb17 = '2blv';
pdbEntries.pdb18 = '2blw';
pdbEntries.pdb19 = '3mnb';
pdbEntries.pdb20 = '3mnc';
pdbEntries.pdb21 = '3mns';
pdbEntries.pdb22 = '3mnx';
pdbEntries.pdb23 = '2j5k';
pdbEntries.pdb24 = '2j5q';
pdbEntries.pdb25 = '2j5r';
%% Loop through each pdb entry and calculate Bdamage for it

%Store all pdb entries in a cell array
pdbCodes = fieldnames(pdbEntries);

%Find the number of pdb entries
numberOfEntries = length(pdbCodes);

%Loop through each pdb entry
for pdb = 1 : numberOfEntries
    %Calculate B damage for each entry
   CalculateBdamage(pdbEntries.(pdbCodes{pdb})) 
end