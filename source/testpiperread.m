function [Energies,nx,ny,nz] = testpiperreadparallel(matfri,parameters)

%%
%% Put the correct size of output

% numcores = 8;
loadlibrary('libMatPiper')
numcores = parameters.piper.numcores;
% Start parallel 
parpool(numcores);


parfor n = 1 : numcores
    % Intialize piper (do only once)
    if ~libisloaded('libMatPiper')
        loadlibrary('libMatPiper')
    end
end
fprintf('Run piper with correct output energy vector \n')

output = int16(0); % Flag true for output of data
%matfri =    0;

nx = parameters.piper.nx;
ny = parameters.piper.ny;
nz = parameters.piper.nz;
piperfile = parameters.piper.piperfileread ;
matlri = matfri + parameters.piper.rotationwidth -1;
if matlri > 69999
    matlri = 69999;
end
FullEnergies = zeros(nx * ny * nz * (matlri - matfri + 1),numcores,'single');
Energies = zeros(x * y * z * (matlri - matfri + 1), 1,'single');

t1 = toc;
[Energies,nx,ny,nz] = calllib('libMatPiper','MatPiper',Energies,1,1,1, matfri, matlri, output, piperfile);
t2 = toc;
fprintf('Total execution time = %f \n', t2 - t1)    
fprintf('\n') 

t1 = toc;

parfor n = 1 : numcores    
    %tic;
    matfri = matfri + (n-1) * parameters.piper.rotationwidth;
    matlri = matfri + parameters.piper.rotationwidth -1;
if matlri > 69999
    matlri = 69999;
end
if matfri > 69999
    matfri = 69999;
end
    Energies = zeros(nx * ny * nz * (matlri - matfri + 1),1,'single');
    fprintf('Start Processor = %d \n', n);    
    [Energies,nx2,ny2,nz2] = calllib('libMatPiper','MatPiper',Energies,1,1,1, matfri, matlri, output, piperfile);
    %toc;
    fprintf('End Processor = %d \n', n);    
    FullEnergies(:,n) = Energies;
    end
%%
t2 = toc;
%unloadlibrary('libMatPiper')

fprintf('Total execution time = %f \n', t2 - t1)    
    fprintf('\n') 



unloadlibrary('libMatPiper')

cd ..


delete(gcp('nocreate'));