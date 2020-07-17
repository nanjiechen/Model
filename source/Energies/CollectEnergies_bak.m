
%Obtain Energies
function [parameters,results] = CollectEnergies_bak(matfri,parameters,results)
l = parameters.stochasticmodel.numofgrids;
n = parameters.piper.nx * parameters.piper.ny * parameters.piper.nz;
d = '/projectnb/uqproj/StochasticDocking/Code/source/piper/input_grids';
files = dir(fullfile(d,'rec_grid*'));
list = {files.name};
m = length(list);
s = ceil(69999/parameters.piper.rotationwidth);
SaveEnergy = parameters.flags.saveenergies;
LoadEnergy = parameters.flags.loadenergies;
%filename = 'energiesMat.mat';
filename = append('energiesMat','_','dim',num2str(parameters.GRFModel.numofeigenvalues),'.mat');
if exist(filename, 'file') == 2 & LoadEnergy == true
    load(filename);
else
% From Rec ( 8 maps x 29)
% For each column of Rec
%   a) Write dx files, it will be 8 files in input_grids/
%   b) Run piper --read-r-grids flag i.e. trypsin_piper_read
%   c) Collect energies (per column)
% We have at the end a energy matrix E (nx ny nz rot) x 29  

EnergiesMat = [];
for j = 1 : l
   for i = 1 : m
   f = parameters.stochasticmodel.RecMap((i-1)*n + 1 : i * n, j);   %acquire each receptor map
   rf = reshape(f,parameters.piper.nx,parameters.piper.ny, parameters.piper.nz); 
   %reshape receptor map to 3D matrix(nx*ny*nz) for dxfile processing
   %name = ["RecMap"+ num2string(i) +".dx"];
   dxfile.densityMatrix = rf;
   dxfile.outfile = list{i};
   cd piper/input_grids
   mat2dx(dxfile);
   cd ..
   cd ..
   end
    [Energies,nx,ny,nz]= testpiperread(matfri,parameters);
  
    EnergiesMat = [EnergiesMat Energies];
end

 if SaveEnergy == true
        save(filename, 'EnergiesMat','-v7.3');
 end
    
end
results.EnergiesMat = EnergiesMat;
end