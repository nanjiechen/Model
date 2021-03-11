

function parameters = piperinitializationsort(parameters)
% Load piper dynamic shared library
% loadlibrary('./piper/libMatPiper')
% Intialize piper (do only once)

% Initialize piper variables
parameters.piper.matfri = 0;
parameters.piper.matlri = 0;
%piper.output = int16(1);
parameters.piper.output = int16(0); % Flag false for do not output the data


parameters.piper.Energies = 1;
parameters.piper.piperfileread = 'piper_trypsin_read'; % Run piper from input dx files
parameters.piper.piperfile     = 'piper_trypsin'; % Run piper from PDB files

parameters.piper.DxFile = '';
parameters.piper.md = '../../data/md_trypsin/';

% Library call to first determine size of the results
cd piper

loadlibrary('libMatPiper')

matfri = 0;
matlri = 0;
output = int16(0); % Flag false for do not output the data

piperfile = 'piper_trypsin';
%piperfile = 'piper_argfile';
%piperfile = 'prot_001_out_nmin.pdb';

Energies = 0;

tic;
fprintf('Calculate size of Cube \n')
t1 = toc;
[Energies,nx,ny,nz] = calllib('libMatPiper','MatPiper',Energies,1,1,1, ...
                                matfri, matlri, output, piperfile);             
t2 = toc;

fprintf('Total execution time = %f \n', t2 - t1)    
fprintf('\n') 

fprintf('Run piper with correct output energy vector: First test \n')

%%
% Intialize sort piper (do only once)
loadlibrary('libMatPiperSort')

output = int16(1); % Flag true for output of data
matfri = 0; % First rotation
%matlri = 0;   % Last rotation
matlri = 10;
Energies = zeros(nx * ny * nz * (matlri - matfri + 1), 1,'single');

nrind = int16(matlri - matfri + 1);
rind = int16(matlri:-1:matfri);

t1 = toc;
[Energies,nx,ny,nz] = calllib('libMatPiperSort','MatPiper',Energies,1,1,1, nrind, rind, output, piperfile);
t2 = toc;
fprintf('Total execution time = %f \n', t2 - t1)    
fprintf('\n')
% Save piper structure in parameters
% parameters.piper = piper;

parameters.piper.nx = nx;
parameters.piper.ny = ny;
parameters.piper.nz = nz;
parameters.piper.numrotations = 10;

% Total memory per single rotation
parameters.piper.numcores = 8;
parameters.piper.memperrot = nx * ny * nz * 4;
parameters.piper.memory = parameters.piper.memory;
%%

cd ..





end