% Run piper to create DX files

cd piper

loadlibrary('libMatPiper')

matfri = 0;
matlri = 0;
output = int16(0); % Flag false for do not output the data

piperfile = parameters.piper.piperfile;
Energies = 1;

% Output size of cube
tic;
fprintf('Calculate size of Cube \n');
t1 = toc;
[Energies,nx,ny,nz] = calllib('libMatPiper','MatPiper',Energies,1,1,1, ...
                                matfri, matlri, output, piperfile);             

t2 = toc;
fprintf('Total execution time = %f \n', t2 - t1)    
fprintf('\n') 

unloadlibrary('libMatPiper')
cd ..
