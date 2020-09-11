
function parameters = piperinitializationDebug(parameters)
% Load piper dynamic shared library
% loadlibrary('./piper/libMatPiper')

% Initialize piper variables
parameters.piper.matfri = 0;
parameters.piper.matlri = 0;
%piper.output = int16(1);
parameters.piper.output = int16(0); % Flag false for do not output the data


parameters.piper.Energies = 1;
parameters.piper.piperfileread = 'piper_trypsin_read_deb'; % Run piper from input dx files
parameters.piper.piperfile     = 'piper_trypsin_read_deb'; % Run piper from PDB files

parameters.piper.DxFile = '';
parameters.piper.md = '../../data/md_trypsin/';

% Library call to first determine size of the results
cd piper

% Intialize piper (do only once)
loadlibrary('libMatPiper')

matfri = 0;
matlri = 0;
output = int16(0); % Flag false for do not output the data

piperfile = 'piper_trypsin_read_deb';

Energies = 1;

tic;
fprintf('Calculate size of Cube \n')
t1 = toc;
[Energies,nx,ny,nz] = calllib('libMatPiper','MatPiper',Energies,1,1,1, ...
                                matfri, matlri, output, piperfile);             
t2 = toc;

fprintf('Total execution time = %f \n', t2 - t1)    
fprintf('\n') 

parameters.piper.nx = nx;
parameters.piper.ny = ny;
parameters.piper.nz = nz;
parameters.piper.numrotations = 10;

% Total memory per single rotation
parameters.piper.numcores = 8;
parameters.piper.memperrot = nx * ny * nz * 4;
parameters.piper.memory = parameters.piper.memory;


% Save piper structure in parameters
% parameters.piper = piper;

%%

cd ..





end