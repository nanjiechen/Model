function [Energies,nx,ny,nz] = testpipersorting(rind,parameters)
% Intialize sort piper (do only once)
cd piper
fprintf('Run piper with correct output energy vector \n')
% loadlibrary('libMatPiperSort')
x = parameters.piper.nx;
y = parameters.piper.ny;
z = parameters.piper.nz;
piperfile = parameters.piper.piperfileread ;
output = int16(1); % Flag true for output of data
nrind = length(rind);
Energies = zeros(x * y * z * nrind, 1,'single');



t1 = toc;
[Energies,nx,ny,nz] = calllib('libMatPiperSort','MatPiper',Energies,1,1,1, nrind, rind, output, piperfile);
t2 = toc;
fprintf('Total execution time = %f \n', t2 - t1)    
fprintf('\n') 
cd ..