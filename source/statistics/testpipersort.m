

%%
%% Put the correct size of output

% Library call to first determine size of the results
cd piper

% Intialize piper (do only once)
loadlibrary('libMatPiperSort')

nrind = int32(1);
rind = int32(1);
output = int32(0); % Flag false for do not output the data

piperfile = 'piper_trypsin';
%piperfile = 'piper_argfile';
%piperfile = 'prot_001_out_nmin.pdb';

Energies = 0;

tic;
fprintf('Calculate size of Cube \n')
t1 = toc;
[Energies,nx,ny,nz] = calllib('libMatPiperSort','MatPiper',Energies,1,1,1, ...
                                nrind, rind, output, piperfile);             
t2 = toc;

fprintf('Total execution time = %f \n', t2 - t1)    
fprintf('\n') 

fprintf('Run piper with correct output energy vector: First test \n')



%%
output = int32(1); % Flag true for output of data
matfri = 1; % First rotation
%matlri = 0;   % Last rotation
matlri = 100;
Energies = zeros(nx * ny * nz * (matlri - matfri + 1), 1,'single');

nrind = int32(matlri - matfri + 1);
rind = int32(matfri:matlri);

t1 = toc;
[Energies,nx,ny,nz] = calllib('libMatPiperSort','MatPiper',Energies,1,1,1, nrind, rind, output, piperfile);
t2 = toc;
fprintf('Total execution time = %f \n', t2 - t1)    
fprintf('\n') 


keyboard


%%
fprintf('Run piper with correct output energy vector: Second test \n')

output = int16(1); % Flag true for output of data
matfri = 0; % First rotation
%matlri = 0;   % Last rotation
matlri2 = 19;
Energies2 = zeros(nx * ny * nz * (matlri2 - matfri + 1), 1,'single');

t1 = toc;
[Energies2,nx,ny,nz] = calllib('libMatPiperSort','MatPiper',Energies2,1,1,1, matfri, matlri2, output, piperfile);
t2 = toc;
fprintf('Total execution time = %f \n', t2 - t1)    
fprintf('\n') 
% en = norm(Energies2(1:length(Energies)) - Energies);
% fprintf('norm of difference = %f \n',en)

fprintf('matlri1 = %f matlri2= %f \n',matlri,matlri2)


unloadlibrary('libMatPiperSort')

cd ..
