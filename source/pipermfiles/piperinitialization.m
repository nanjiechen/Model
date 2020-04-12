
function parameters = piperinitialization(parameters)
% Load piper dynamic shared library
% loadlibrary('./piper/libMatPiper')

% Initialize piper variables
piper.matfri = 0;
piper.matlri = 0;
%piper.output = int16(1);
piper.output = int16(0); % Flag false for do not output the data


piper.Energies = 1;
piper.piperfileread = 'piper_trypsin_read'; % Run piper from input dx files
piper.piperfile     = 'piper_trypsin'; % Run piper from PDB files

piper.DxFile = '';
piper.md = '../../data/md_trypsin/';

% Library call to first determine size of the results
cd piper

% Intialize piper (do only once)
loadlibrary('libMatPiper')

matfri = 0;
matlri = 0;
output = int16(0); % Flag false for do not output the data

piperfile = 'piper_trypsin';

Energies = 1;

tic;
fprintf('Calculate size of Cube \n')
t1 = toc;
[Energies,nx,ny,nz] = calllib('libMatPiper','MatPiper',Energies,1,1,1, ...
                                matfri, matlri, output, piperfile);             
t2 = toc;

fprintf('Total execution time = %f \n', t2 - t1)    
fprintf('\n') 

piper.nx = nx;
piper.ny = ny;
piper.nz = nz;
piper.numrotations = 10;

% Total memory per single rotation
piper.memperrot = nx * ny * nz * 4;
piper.memory = parameters.piper.memory;

% Save piper structure in parameters
parameters.piper = piper;

%%

cd ..





end
