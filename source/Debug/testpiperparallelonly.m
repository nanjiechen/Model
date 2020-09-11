

%%
%% Put the correct size of output

% Library call to first determine size of the results
cd piper

numcores = 6;
loadlibrary('libMatPiper')

% Start parallel 
parpool(numcores);


parfor n = 1 : numcores
    % Intialize piper (do only once)
    if ~libisloaded('libMatPiper')
        loadlibrary('libMatPiper')
    end
end

matfri = 0;
matlri = 0;
output = int16(0); % Flag false for do not output the data

piperfile = 'piper_argfile';

Energies = 1;

tic;
fprintf('Calculate size of Cube \n')
t1 = toc;
[Energies,nx,ny,nz] = calllib('libMatPiper','MatPiper',Energies,1,1,1, ...
                                matfri, matlri, output, piperfile);             
t2 = toc;

fprintf('Total execution time = %f \n', t2 - t1)    
fprintf('\n') 

fprintf('Run piper with correct output energy vector \n')

output = int16(1); % Flag true for output of data
matfri = 0; % First rotation
%matlri = 0;   % Last rotation
matlri = 100;


%%
FullEnergies = zeros(nx * ny * nz * (matlri - matfri + 1),numcores,'single');
%nx2 = zeros(numcores,1);
%ny2 = nx2;
%nz2 = nx2;


Energies = zeros(nx * ny * nz * (matlri - matfri + 1), 1,'single');
t1 = toc;
[Energies,nx,ny,nz] = calllib('libMatPiper','MatPiper',Energies,1,1,1, matfri, matlri, output, piperfile);
t2 = toc;
fprintf('Total execution time = %f \n', t2 - t1)    
fprintf('\n') 


t1 = toc;

parfor n = 1 : numcores    
    %tic;
    Energies = zeros(nx * ny * nz * (matlri - matfri + 1),1,'single');
    fprintf('Start Processor = %d \n', n);    
    [Energies,nx2,ny2,nz2] = calllib('libMatPiper','MatPiper',Energies,1,1,1, matfri, matlri, output, piperfile);
    %toc;
    fprintf('End Processor = %d \n', n);    
    FullEnergies(:,n) = Energies;
    end
%%
t2 = toc;

fprintf('Total execution time = %f \n', t2 - t1)    
    fprintf('\n') 



%unloadlibrary('libMatPiper')

cd ..


delete(gcp('nocreate'));