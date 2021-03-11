function [Energies,nx,ny,nz] = testpiperreadbak(matfri,parameters)

%%
%% Put the correct size of output

cd piper
fprintf('Run piper with correct output energy vector \n')

output = int16(0); % Flag true for output of data
%matfri =    0;

x = parameters.piper.nx;
y = parameters.piper.ny;
z = parameters.piper.nz;
piperfile = parameters.piper.piperfileread ;
matlri = matfri + parameters.piper.rotationwidth -1;
if matlri > 69999
    matlri = 69999;
end

Energies = zeros(x * y * z * (matlri - matfri + 1), 1,'single');

t1 = toc;
[Energies,nx,ny,nz] = calllib('libMatPiper','MatPiper',Energies,1,1,1, matfri, matlri, output, piperfile);
t2 = toc;
fprintf('Total execution time = %f \n', t2 - t1)    
fprintf('\n') 

%unloadlibrary('libMatPiper')

cd ..