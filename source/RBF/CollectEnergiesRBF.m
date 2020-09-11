function [parameters,results] = CollectEnergiesRBF(parameters,results,methods)


s = size(results.SortedMat,1);
z = parameters.piper.nx * parameters.piper.ny * parameters.piper.nz;
rotations = results.SortedMat(:,1) -1;
parameters.piper.rotationwidth  = 1;
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
Trans = results.SortedMat(:,2);

for k = 1:s
      matfri = rotations(k);
      ind = Trans(k);
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
    en = Energies(ind);
    Block = [Block en];
       clear Energies;
     save('../data/PiperData/xyzDataMat','-v7.3');

end
results.Block = Block;
  [parameters,results] = methods.RBF(parameters,results);
  ys = results.RBFys';
  RBFapprox = [RBFapprox ; ys];
       EnergiesMat = [EnergiesMat;Block]; 
       
 end
end


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