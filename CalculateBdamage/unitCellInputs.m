%% Function to create the inputs necessary to run PDBCUR
function inputs = unitCellInputs()

%write the name of the input file that will be used in PDBCUR run.
inputs.filename = 'pdbcurInputUnitCellGeneration.txt';
%input to remove hydrogen atoms from file
inputs.delhydrogen = 'delhydrogen';
%input to keep the most probable alternate conformation
inputs.mostprob = 'mostprob';
%input to remove anisotropic (ANISOU) records from the file
inputs.noanisou = 'noanisou';
%input to generate the unit cell from the symmetry operations
inputs.genunit = 'genunit';

end

