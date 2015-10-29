%% Function which calls the program PDBCUR to process the pdb file
function processPDB(inputPDB,outputPDB,inputParameters)

%Create name for the input file
inputFileNameForPBDCUR = inputParameters.filename;

%Create the input file for PDBCUR
createInputFileForPDBCUR(inputFileNameForPBDCUR,inputParameters); 

%Create the system command to run PDBCUR
runPDBCURCommand = sprintf('pdbcur xyzin %s xyzout %s < %s',inputPDB,outputPDB,inputFileNameForPBDCUR);

%Run PDBCUR to process the pdb file
system(runPDBCURCommand,'-echo');

end

