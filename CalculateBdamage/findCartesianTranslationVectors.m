%% Function to convert from unit cell parameters to cartesian coordinates
function [aCartesianVector,bCartesianVector,cCartesianVector] = findCartesianTranslationVectors(unitCellParams)

%Explicitly write the params of the Unit Cell
a       = unitCellParams(1);
b       = unitCellParams(2);
c       = unitCellParams(3);

%convert the angles from degrees to radians
alpha   = degtorad(unitCellParams(4));
beta    = degtorad(unitCellParams(5));
gamma   = degtorad(unitCellParams(6));

%Define parameter v which is the final input of the conversion matrix below
v = sqrt(1 - (cos(alpha))^2 - (cos(beta))^2 - (cos(gamma))^2 ...
    + 2 * cos(alpha) * cos(beta) * cos(gamma));

%Define the elements of the conversion matrix
a11 = a;
a12 = b * cos(gamma);
a13 = c * cos(beta);

a21 = 0;
a22 = b * sin(gamma);
a23 = c * ((cos(alpha) - (cos(beta) * cos(gamma)))/(sin(gamma)));

a31 = 0;
a32 = 0;
a33 = c * (v / sin(gamma));

%Create the Matrix
conversionMatrix = [a11, a12, a13; ...
    a21, a22, a23; ...
    a31, a32, a33];

%define the fractional basis vectors for each direction a, b and c
aBasisVector = [1;0;0];
bBasisVector = [0;1;0];
cBasisVector = [0;0;1];

%Convert the lattice basis vector in each direction to Cartesian
%coordinates
aCartesianVector = conversionMatrix * aBasisVector;
bCartesianVector = conversionMatrix * bBasisVector;
cCartesianVector = conversionMatrix * cBasisVector;

end