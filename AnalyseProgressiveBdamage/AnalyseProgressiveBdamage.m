%% Script to analyse how Bdamage varies in a progressive data set
%Input: file path to folder containing original and subsequent data sets 
%(in .txt files) for which Bdamage values have already been calculated
%Copyright 2015 Thomas Dixon
fprintf('Copyright 2015 Thomas Dixon\n')
fprintf('\n')
function AnalyseProgressiveBdamage(pathToBdamageFolder)
%% User Input selection of tasks to perform

%choose whether to present Bdamage values on a single scatter graph for all
%runs
presentBd = true;

%choose whether to calculate and output Bprogression
calcBprogression = ~true;

%% Inform the user of the time of initialisation of the script
%Start timer to calculate time taken to run the script
mainTimer = tic;

fprintf('################################################################\n')
fprintf('################################################################\n')
fprintf('################################################################\n')
fprintf('################## Program to analyse B Damage #################\n')
fprintf('################################################################\n')
fprintf('\n\n')

%Print the time and date that the program was started
%Get the time and date for the program 
programStartTime = clock;
theYear = programStartTime(1);
theMonth = programStartTime(2);
theDay = programStartTime(3);
theHour = programStartTime(4);
theMinute = programStartTime(5);
theSecond = programStartTime(6);

fprintf('This program was run on %d/%d/%d at %02d:%02d:%02.0f.\n\n',theDay,theMonth,theYear,theHour,theMinute,theSecond);

%% Extract Input File Series
fprintf('****************************************************************\n')
fprintf('************************ Input Section *************************\n')
fprintf('\n')

%Create a struct containg the filenames of the files within the input
%folder that end in .txt
fileStruct = dir(fullfile(sprintf('%s',pathToBdamageFolder),'*.txt'));

%Copy the Bdamage files to the workspace
numberOfFiles = length(fileStruct);
for file = 1 : numberOfFiles
    filename = fileStruct(file).name;
    copyfile(sprintf('%s\\%s',pathToBdamageFolder,filename),sprintf('set%02d.txt',file));
    fprintf(sprintf('set%02d.txt created\n',file))
end

fprintf('\n')
fprintf('******************* End of Input Section ***********************\n')
fprintf('****************************************************************\n')
fprintf('\n')
fprintf('----------------------------------------------------------------\n')
fprintf('\n')

%% Loop through each Bdamage file and extract Info
fprintf('****************************************************************\n')
fprintf('**************** Parsing all Files Section *********************\n')
fprintf('\n')

%Create a Cell Array the size of the number of files
allBdamage = cell(1,numberOfFiles);

%Find the number of atoms in the file

%Read in the Bdamage file
file = 1;
bdamageFile = fileread(sprintf('set%02d.txt',file));

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

%Loop through each Bdamage File and parse the Bdamage file
for file = 1 : numberOfFiles
    fprintf('dataArray%02d\n',file);
    allBdamage{file} = parseBdamage(sprintf('set%02d.txt',file)); 
end

fprintf('\n')
fprintf('****************** End of Parsing Section **********************\n')
fprintf('****************************************************************\n')
fprintf('\n')
fprintf('----------------------------------------------------------------\n')
fprintf('\n')

%% Present Bdamage Values for each Data Set
fprintf('****************************************************************\n')
fprintf('**************** Data Presentation Section *********************\n')
fprintf('\n')

%Only if presentBdamage set to true
if presentBd
    presentBdamage(allBdamage,numberOfFiles,pathToBdamageFolder)  
end

%If user selected not to present the data, output message
if ~presentBd
    fprintf('User selected not to present Bdamage')
    fprintf('\n')
end

fprintf('\n')
fprintf('************** End of Data Presentation Section ****************\n')
fprintf('****************************************************************\n')
fprintf('\n')
fprintf('----------------------------------------------------------------\n')
fprintf('\n')
%% Calculate Bprogression
fprintf('****************************************************************\n')
fprintf('****************** Calculate Bprogression **********************\n')
fprintf('\n')

%Only if calcBprogression set to true
if calcBprogression
    [allBdamage] = calculateBprogression(allBdamage,numberOfFiles,numberOfAtoms,pathToBdamageFolder);
end

%If user selected not to claculate Bprogression, output message
if ~calcBprogression
    fprintf('User selected not to calculate Bprogression')
    fprintf('\n')
end

fprintf('\n')
fprintf('**************** End of calculate Bprogression *****************\n')
fprintf('****************************************************************\n')
fprintf('\n')
fprintf('----------------------------------------------------------------\n')
fprintf('\n')

% %% Perform Numeric Analysis of Bdamage Values
% fprintf('****************************************************************\n')
% fprintf('**************** Data Presentation Section *********************\n')
% fprintf('\n')

%% Delete the temporary 'setxx' files

for file = 1 : numberOfFiles
    delete(sprintf('set%02d.txt',file));
end

%% Stop timing the time taken for the program to run
timeTaken = toc(mainTimer);

fprintf('\n');
fprintf('Total time taken for script to run was %.0f minutes and %.0f seconds.\n',floor(timeTaken/60),rem(timeTaken,60));

end