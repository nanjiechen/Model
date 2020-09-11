function [parameters,results] = CollectEnergiesxyz(parameters,results)

% List of pdb files
s = ceil(69999/parameters.piper.rotationwidth);
z = parameters.piper.nx * parameters.piper.ny * parameters.piper.nz;

LoadEnergy = parameters.flags.loadenergies;
SaveEnergy = parameters.flags.saveenergies;
 filename = append('EnergiesMat','_',num2str(parameters.GRFModel.numofeigenvalues),'xyz','.mat');
% filename = 'EnergiesMat_2_level4_0.1.mat';
if exist(filename, 'file') == 2 & LoadEnergy == true
    load(filename);
else
EnergiesMat = [];
matfri = 0;
d = '/projectnb/uqproj/StochasticDocking/Code/data/PDBApprox/';
%files = dir(['/projectnb/uqproj/StochasticDocking/Code/data/md_trypsin/','*.pdb']);
f = dir([d,'*approximation*.pdb']);
%files =dir([parameters.piper.md,'/*.pdb']);

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

%  for k = 942:942
for k = 1:s
      matfri = (k -1) * parameters.piper.rotationwidth;
     if matfri + parameters.piper.rotationwidth > 70000
    BlockIndex = results.TransIndex( (k - 1) * parameters.piper.rotationwidth + 1 : end);  
    else
    BlockIndex = results.TransIndex( (k - 1) * parameters.piper.rotationwidth + 1 : parameters.piper.rotationwidth * k);
    end
    K = BlockIndex - z * (k -1) * parameters.piper.rotationwidth;
      Block = [];

for j = 1:n
   % p = append(parameters.piperdf -H.md,list{j});
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
     save('../data/PiperData/xyzDataMat','-v7.3');

end
    matfri = k * parameters.piper.rotationwidth;
       EnergiesMat = [EnergiesMat;Block]; 
       
 end
end


[tf, Index_order]= ismember(results.OrigRotation,results.SortedRotation);
% [results.OrigRotation, Index_order] = sort(results.SortedRotation);
EnergiesMat = EnergiesMat(Index_order,:);

fname = append('../data/PiperData/','EnergiesMat','_',num2str(parameters.GRFModel.numofeigenvalues),'_xyz_ordered');
if SaveEnergy == true
save (fname, 'EnergiesMat','-v7.3');
end
results.EnergiesMat = EnergiesMat;
 allMin = min(results.EnergiesMat, [],'all');
  results.EnergyMin = allMin;
 if results.EnergyMin < 1
    results.ScaledEnergies = results.EnergiesMat + ( 1 - results.EnergyMin);
    else
    results.ScaledEnergies = results.EnergiesMat;
 end
%  cd piper
%  [status,cmdout] = system('> multirec_trypsin_expand.pdb');
%  cd ..
end
    