function [dataArray] = parseBdamage(BdamageFilePath)

%% Read in the Bdamage file

bdamageFile = fileread(BdamageFilePath);

fprintf('Now extracting Atom and Bdamage information.\n')

%% Parse the file to get table containing HETATM info

%Create an expression in the pdb file that starts with the word "ATOM"
%then any number of white spaces then a number then any number of white
%spaces and then a capital letter between A and Z, then any number of white
%spaces then a capital letter between A and Z.
expressionToFind = 'ATOM\s*\d\s*';

%Find the expression in the pdb file and return the indices where it has
%been found.
startIndexOfATOMbdamageInfo = regexp(bdamageFile,expressionToFind);

%Calculate the number of atoms in the file
numberOfAtoms = length(startIndexOfATOMbdamageInfo);

%Create an expression in the pdb file that starts with the word "HETATM"
%then any number of white spaces then a number then any number of white
%spaces and then a capital letter between A and Z, then any number of white
%spaces then a capital letter between A and Z.
expressionToFind = 'HETATM\s*\d\s*';

%Find the expression in the pdb file and return the indices where it has
%been found.
startIndexOfHETATMbdamageInfo = regexp(bdamageFile,expressionToFind);

%Calculate the number of Hetero atoms in the file
numberOfHetatms = length(startIndexOfHETATMbdamageInfo);

%% Create a cell array where each row contains information about each atom

%preallocate cell array to store the information
dataArray = cell(numberOfAtoms + numberOfHetatms,21);

%Create a single string containing all of the rows as a table by
%concatenating the strings in the loop for each atom
tableOfInfo = '';
for eachRow = 1:length(startIndexOfATOMbdamageInfo)
    atomBdamInfo = bdamageFile(startIndexOfATOMbdamageInfo(eachRow) : startIndexOfATOMbdamageInfo(eachRow)+123);
    dataArray{eachRow,1} = atomBdamInfo(1:6);
    dataArray{eachRow,2} = atomBdamInfo(7:11);
    dataArray{eachRow,3} = atomBdamInfo(13:16);
    dataArray{eachRow,4} = atomBdamInfo(17);
    dataArray{eachRow,5} = atomBdamInfo(18:20);
    dataArray{eachRow,6} = atomBdamInfo(22);
    dataArray{eachRow,7} = atomBdamInfo(23:26);
    dataArray{eachRow,8} = atomBdamInfo(27);
    dataArray{eachRow,9} = atomBdamInfo(31:38);
    dataArray{eachRow,10} = atomBdamInfo(39:46);
    dataArray{eachRow,11} = atomBdamInfo(47:54);
    dataArray{eachRow,12} = atomBdamInfo(55:60);
    dataArray{eachRow,13} = atomBdamInfo(61:66);
    dataArray{eachRow,14} = atomBdamInfo(77:78);
    %dataArray{eachRow,15} = atomBdamInfo(79:80);
    dataArray{eachRow,16} = atomBdamInfo(79:84);
    dataArray{eachRow,17} = atomBdamInfo(86:92);
    dataArray{eachRow,18} = atomBdamInfo(94:97);
    dataArray{eachRow,19} = atomBdamInfo(98:100);
    dataArray{eachRow,20} = atomBdamInfo(106:120);
    dataArray{eachRow,21} = atomBdamInfo(122:124);
    tableOfInfo = sprintf('%s%s\n',tableOfInfo,atomBdamInfo(1:end-1));
end

%Now add the Hetero Atoms on the end
% for eachRow = 1:length(startIndexOfHETATMbdamageInfo)
%     atomBdamInfo = bdamageFile(startIndexOfHETATMbdamageInfo(eachRow) : startIndexOfHETATMbdamageInfo(eachRow)+123);
%     dataArray{eachRow + numberOfAtoms,1} = atomBdamInfo(1:6);
%     dataArray{eachRow + numberOfAtoms,2} = atomBdamInfo(7:11);
%     dataArray{eachRow + numberOfAtoms,3} = atomBdamInfo(13:16);
%     dataArray{eachRow + numberOfAtoms,4} = atomBdamInfo(17);
%     dataArray{eachRow + numberOfAtoms,5} = atomBdamInfo(18:20);
%     dataArray{eachRow + numberOfAtoms,6} = atomBdamInfo(22);
%     dataArray{eachRow + numberOfAtoms,7} = atomBdamInfo(23:26);
%     dataArray{eachRow + numberOfAtoms,8} = atomBdamInfo(27);
%     dataArray{eachRow + numberOfAtoms,9} = atomBdamInfo(31:38);
%     dataArray{eachRow + numberOfAtoms,10} = atomBdamInfo(39:46);
%     dataArray{eachRow + numberOfAtoms,11} = atomBdamInfo(47:54);
%     dataArray{eachRow + numberOfAtoms,12} = atomBdamInfo(55:60);
%     dataArray{eachRow + numberOfAtoms,13} = atomBdamInfo(61:66);
%     dataArray{eachRow + numberOfAtoms,14} = atomBdamInfo(77:78);
%     %dataArray{eachRow + numberOfAtoms,15} = atomBdamInfo(79:80);
%     dataArray{eachRow + numberOfAtoms,16} = atomBdamInfo(79:84);
%     dataArray{eachRow + numberOfAtoms,17} = atomBdamInfo(86:92);
%     dataArray{eachRow + numberOfAtoms,18} = atomBdamInfo(93:96);
%     dataArray{eachRow + numberOfAtoms,19} = atomBdamInfo(98:100);
%     dataArray{eachRow + numberOfAtoms,20} = atomBdamInfo(106:120);
%     dataArray{eachRow + numberOfAtoms,21} = atomBdamInfo(122:124);
%     tableOfInfo = sprintf('%s%s\n',tableOfInfo,atomBdamInfo(1:end-1));
% end

fprintf('Finished reading file.\n')

end