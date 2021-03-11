
function [parameters,results] = CollectEnergiesNewTechSort(parameters,results,methods)
parameters.piper.rotationwidth = 1000;
s = ceil(69999/parameters.piper.rotationwidth);
% s = length(results.TransIndex);
z = parameters.piper.nx * parameters.piper.ny * parameters.piper.nz;

LoadEnergy = parameters.flags.loadenergies;
SaveEnergy = parameters.flags.saveenergies;
 filename = append('EnergiesMat','_',num2str(parameters.GRFModel.numofeigenvalues),'xyz','.mat');

if exist(filename, 'file') == 2 & LoadEnergy == true
    load(filename);
else
EnergiesMat = [];
RBFapprox = [];
matfri = 0;
d = '/projectnb/uqproj/StochasticDocking/Code/data/PDBApprox/';
f = dir([d,'*approximation*.pdb']);

list = {f.name};

n = parameters.stochasticmodel.numofgrids;
md = '../../data/PDBApprox/';

fid = fopen(parameters.piper.piperfileread,'r');
i = 1;
tline = fgetl(fid);
A{i}= tline;
 while ischar(tline)
     i = i+1;
    tline = fgetl(fid);
    A{i} = tline;
  
   end
fclose(fid);
EnergiesMat = [];

transitions = results.SortedMat(1:parameters.piper.rotationwidth,2);
rind = results.SortedMat(1:parameters.piper.rotationwidth,1) -1;
% transitions = results.AscendingSortedMat(1:parameters.piper.rotationwidth,2);
% rind = [1:parameters.piper.rotationwidth] -1;
nrind = length(rind);
K = transitions + (z*([1:nrind]-1))';
RBFMat = [];

Block = [];

for j = 1:n

   p = [md, list{j}];
    A{19} = p;

   
    fid = fopen(parameters.piper.piperfileread, 'r+');
    
    for I = 1:22
        fprintf(fid, '%s\n', A{I});
    end
    fclose(fid);
   
   
       
     
    [Energies,nx,ny,nz]= testpipersorting(rind,parameters);
     BlockEnergies = Energies(K);
     Block = [Block BlockEnergies]; 
   
       clear Energies;
%        save('../data/PiperData/RBFDataMatSort','-v7.3');  
end
results.Block = Block;
EnergiesMat = [EnergiesMat;Block]; 
save('../data/PiperData/RBFDataMatSort2','-v7.3');       

end
save('../data/PiperData/NewTechEngSort','-v7.3')


   


% 
% [tf, Index_order]= ismember(results.OrigRotation,results.SortedRotation);
% EnergiesMat = EnergiesMat(Index_order,:);

fname = append('../data/PiperData/','EnergiesMat','_',num2str(parameters.GRFModel.numofeigenvalues),'_',num2str(parameters.GRFModel.KLRescale));if SaveEnergy == true
save (fname, 'EnergiesMat','-v7.3');
end
results.RBFapprox = RBFapprox;
 allMin = min(results.RBFapprox, [],'all');
  results.EnergyMin = allMin;
 if results.EnergyMin < 1
    results.ScaledEnergies = results.RBFapprox + ( 1 - allMin);
    else
    results.ScaledEnergies = results.RBFapprox;
 end

end