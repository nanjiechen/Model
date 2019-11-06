
% Load piper dynamic library
cd piper
loadlibrary('libMatPiper')

% Library call to first determine size of the results
matfri = 0;
matlri = 0;
output = int16(0); % Flag false for do not output the data
%piperfile = 'piper_argfile';
piperfile = 'piper_trypsin';
Energies = 1;

% Output size of cube
tic;
fprintf('Running piper \n')
t1 = toc;
[Energies,nx,ny,nz] = calllib('libMatPiper','MatPiper',Energies,1,1,1, ...
                                matfri, matlri, output, piperfile);             

t2 = toc;
fprintf('Total execution time = %f \n', t2 - t1)    
fprintf('\n') 

% Unload piper library
unloadlibrary('libMatPiper')

cd ..
