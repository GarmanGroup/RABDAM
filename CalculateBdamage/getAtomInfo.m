function [dataArray, unitCellInfo, pdbPreamble, pdbEOF] = getAtomInfo(pdbFilePath)

%% Read in the pdb file

pdbFile = fileread(pdbFilePath);

fprintf('Finished reading file.\n')
fprintf('Now extracting Unit Cell Information and Atomic coordinates.\n')

%% Find the Unit Cell dimension

%Create expression to find the line in the PDB file that corresponds to the
%Unit Cell data
expressionToFind = 'CRYST1';

%Find the Unit Cell information in the pdb file
indexForUnitCellInfo = regexp(pdbFile,expressionToFind);

%Split the string containing the unit cell information
UnitCellInfoAsText = strsplit(pdbFile(indexForUnitCellInfo : indexForUnitCellInfo + 80));

%Preallocate vector to store unit cell information
unitCellInfo = zeros(1,6);

%Store unit cell information in a vector. The information stored is:
%Unit Cell dimension: a
%Unit Cell dimension: b
%Unit Cell dimension: c
%Unit Cell angle: alpha
%Unit Cell angle: beta
%Unit Cell angle: gamma
unitCellInfo(1) = str2double(UnitCellInfoAsText{2});
unitCellInfo(2) = str2double(UnitCellInfoAsText{3});
unitCellInfo(3) = str2double(UnitCellInfoAsText{4});
unitCellInfo(4) = str2double(UnitCellInfoAsText{5});
unitCellInfo(5) = str2double(UnitCellInfoAsText{6});
unitCellInfo(6) = str2double(UnitCellInfoAsText{7});

%% Parse the File to get table containing ATOM and HETATM info

%Create an expression in the pdb file that starts with the word "ATOM"
%then any number of white spaces then a number then any number of white
%spaces and then a capital letter between A and Z, then any number of white
%spaces then a capital letter between A and Z.
expressionToFind = 'ATOM\s*\d\s*';

%Find the expression in the pdb file and return the indices where it has
%been found.
startIndexOfAtomInfo = regexp(pdbFile,expressionToFind);

%Create an expression in the pdb file that starts with the word "HETATM"
%then any number of white spaces then a number then any number of white
%spaces and then a capital letter between A and Z, then any number of white
%spaces then a capital letter between A and Z.
expressionToFind = 'HETATM\s*\d\s*';

%Find the expression in the pdb file and return the indices where it has
%been found.
startIndexOfHETATMInfo = regexp(pdbFile,expressionToFind);

%% Create a cell array where each row contains information about each atom

%preallocate cell array to store the information
dataArray = cell(length(startIndexOfAtomInfo) + length(startIndexOfHETATMInfo),15);

%Create a single string containing all of the rows as a table by
%concatenating the strings in the loop for each atom
for eachRow = 1:length(startIndexOfAtomInfo)
    atomInfo = pdbFile(startIndexOfAtomInfo(eachRow) : startIndexOfAtomInfo(eachRow)+78);
    dataArray{eachRow,1} = atomInfo(1:6);
    dataArray{eachRow,2} = atomInfo(7:11);
    dataArray{eachRow,3} = atomInfo(13:16);
    dataArray{eachRow,4} = atomInfo(17);
    dataArray{eachRow,5} = atomInfo(18:20);
    dataArray{eachRow,6} = atomInfo(22);
    dataArray{eachRow,7} = atomInfo(23:26);
    dataArray{eachRow,8} = atomInfo(27);
    dataArray{eachRow,9} = atomInfo(31:38);
    dataArray{eachRow,10} = atomInfo(39:46);
    dataArray{eachRow,11} = atomInfo(47:54);
    dataArray{eachRow,12} = atomInfo(55:60);
    dataArray{eachRow,13} = atomInfo(61:66);
    dataArray{eachRow,14} = atomInfo(77:78);
    %dataArray{eachRow,15} = atomInfo(79:80);
end

%Now add the Heavy Atoms on the end
for eachRow = 1:length(startIndexOfHETATMInfo)
    atomInfo = pdbFile(startIndexOfHETATMInfo(eachRow) : startIndexOfHETATMInfo(eachRow)+78);
    dataArray{eachRow + length(startIndexOfAtomInfo),1} = atomInfo(1:6);
    dataArray{eachRow + length(startIndexOfAtomInfo),2} = atomInfo(7:11);
    dataArray{eachRow + length(startIndexOfAtomInfo),3} = atomInfo(13:16);
    dataArray{eachRow + length(startIndexOfAtomInfo),4} = atomInfo(17);
    dataArray{eachRow + length(startIndexOfAtomInfo),5} = atomInfo(18:20);
    dataArray{eachRow + length(startIndexOfAtomInfo),6} = atomInfo(22);
    dataArray{eachRow + length(startIndexOfAtomInfo),7} = atomInfo(23:26);
    dataArray{eachRow + length(startIndexOfAtomInfo),8} = atomInfo(27);
    dataArray{eachRow + length(startIndexOfAtomInfo),9} = atomInfo(31:38);
    dataArray{eachRow + length(startIndexOfAtomInfo),10} = atomInfo(39:46);
    dataArray{eachRow + length(startIndexOfAtomInfo),11} = atomInfo(47:54);
    dataArray{eachRow + length(startIndexOfAtomInfo),12} = atomInfo(55:60);
    dataArray{eachRow + length(startIndexOfAtomInfo),13} = atomInfo(61:66);
    dataArray{eachRow + length(startIndexOfAtomInfo),14} = atomInfo(77:78);
    %dataArray{eachRow + length(startIndexOfAtomInfo),15} = atomInfo(79:80);
end

%% Get start and end of file 
%Get the string that corresponds to the entire PDB file up until the ATOM
%records
pdbPreamble = pdbFile(1:startIndexOfAtomInfo-1);

%Get the string that corresponds to the end the pdb file from the end of
%the atomic coordinates
pdbEOF = pdbFile(startIndexOfHETATMInfo(eachRow)+81 : end);
end

