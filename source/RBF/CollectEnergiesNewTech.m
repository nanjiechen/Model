
function [parameters,results] = CollectEnergiesNewTech(parameters,results,methods)
parameters.piper.rotationwidth = 3000;
s = ceil(69999/parameters.piper.rotationwidth);
% s = length(results.TransIndex);

z = parameters.piper.nx * parameters.piper.ny * parameters.piper.nz;
rotations = results.SortedMat(:,1) -1;



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
results.RBFCoeffVec = [];
Trans = results.SortedMat(:,2);
RBFMat = [];
for k = 1:s
%       matfri = rotations(k);
%       ind = Trans(k);
  matfri = (k -1) * parameters.piper.rotationwidth;
     if matfri + parameters.piper.rotationwidth > 70000
    BlockIndex = results.TransIndex( (k - 1) * parameters.piper.rotationwidth + 1 : end);  
    else
    BlockIndex = results.TransIndex( (k - 1) * parameters.piper.rotationwidth + 1 : parameters.piper.rotationwidth * k);
    end
    K = BlockIndex - z * (k -1) * parameters.piper.rotationwidth;
      Block = [];

for j = 1:n

   p = [md, list{j}];
    A{19} = p;

   
    fid = fopen(parameters.piper.piperfileread, 'r+');
    
    for I = 1:22
        fprintf(fid, '%s\n', A{I});
    end
    fclose(fid);
   
   
       
     
    [Energies,nx,ny,nz]= testpiperread(matfri,parameters);
     BlockEnergies = Energies(K);
     Block = [Block BlockEnergies]; 
       clear Energies;
     save('../data/PiperData/RBFDataMat','-v7.3');
     
%     en = Energies(ind);
%     Block = [Block en];
%        clear Energies;
%      save('../data/PiperData/xyzDataMat','-v7.3');

% results.Block = Block;
%   [parameters,results] = methods.RBF(parameters,results);
%   ys = results.RBFys';
%   RBFapprox = [RBFapprox ; ys];
end
results.Block = Block;
EnergiesMat = [EnergiesMat;Block]; 
matfri = k * parameters.piper.rotationwidth;
       
end
 save('../data/PiperData/NewTechEng','-v7.3')
end
 load('Weihgts_roots.mat'); % the data file name may vary!
         
         % Notice that you should already constructed an object before
             R.roots_1d=roots;
             R.weights_1d=weights;
             quabak=R.Baking;

   



% [tf, Index_order]= ismember(results.OrigRotation,results.SortedRotation);
% EnergiesMat = EnergiesMat(Index_order,:);

fname = append('../data/PiperData/','EnergiesMat','_',num2str(parameters.GRFModel.numofeigenvalues),'_RBF_ordered');
if SaveEnergy == true
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