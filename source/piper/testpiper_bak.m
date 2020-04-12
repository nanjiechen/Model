

%%
%% Put the correct size of output

cd piper
fprintf('Run piper with correct output energy vector \n')

output = int16(1); % Flag true for output of data
matfri =    5445;

nx = parameters.piper.nx;
ny = parameters.piper.ny;
nz = parameters.piper.nz;
piperfile = parameters.piper.piperfile;
matlri = matfri + parameters.piper.numrotations;

Energies = zeros(nx * ny * nz * (matlri - matfri + 1), 1,'single');

t1 = toc;
[Energies,nx,ny,nz] = calllib('libMatPiper','MatPiper',Energies,1,1,1, matfri, matlri, output, piperfile);
t2 = toc;
fprintf('Total execution time = %f \n', t2 - t1)    
fprintf('\n') 

%unloadlibrary('libMatPiper')

cd ..
