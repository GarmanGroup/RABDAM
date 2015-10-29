%% Script to calculate Bdamage values for multpile input pdbs contained in a single folder
%Input: file path to folder contiaining all data sets as .pdb files
function inputProgressiveDataset(pathToFolder)

%% Create Input file series

%Create a struct containg the filenames of the files within the input
%folder that end in .pdb
fileStruct = dir(fullfile(sprintf('%s',pathToFolder),'*.pdb'));

%Create a puppet filepath for each of the files
numberOfFiles = length(fileStruct);
for file = 1 : numberOfFiles
    filename = fileStruct(file).name;
    filepath = sprintf('%s\\%s',pathToFolder,filename);
    
    %Calculate Bdamage for the file referenced by this filepath
    CalculateBdamage(sprintf('%s',filepath))
end

end