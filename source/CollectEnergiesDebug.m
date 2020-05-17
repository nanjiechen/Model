
%Obtain Energies
function [parameters,results] = CollectEnergiesDebug(matfri,results,parameters)
parameters.piper.rotationwidth = floor(parameters.piper.memory / ...
    (parameters.piper.memperrot * 100));
matfri = 36222;
l = parameters.stochasticmodel.numofgrids;
n = parameters.piper.nx * parameters.piper.ny * parameters.piper.nz;
inputd = '/projectnb/uqproj/StochasticDocking/Code/source/piper/input_grids';
files = dir(fullfile(inputd,'rec_grid*'));
list = {files.name};
m = length(list);
EnergiesMat = [];

s = ceil(69999/parameters.piper.rotationwidth);
% List of pdb files
data = [];
d = '/projectnb/uqproj/StochasticDocking/Code/data/md_trypsin/';
pdbfiles = dir([d,'*out*.pdb']);
pdbls = {pdbfiles.name};
pdblength = length(pdbfiles);

SaveEnergy = parameters.flags.saveenergies;
LoadEnergy = parameters.flags.loadenergies;
B = results.spatial2ind;
%filename = 'energiesMat.mat';
filename = append('energiesMat','_','dim',num2str(parameters.GRFModel.numofeigenvalues),'_deb','.mat');
if exist(filename, 'file') == 2 & LoadEnergy == true
    load(filename);
else
% From Rec ( 8 maps x 29)
% For each column of Rec
%   a) Write dx files, it will be 8 files in input_grids/
%   b) Run piper --read-r-grids flag i.e. trypsin_piper_read
%   c) Collect energies (per column)
% We have at the end a energy matrix E (nx ny nz rot) x 29  


fid = fopen('piper_trypsin_read_deb','r');
i = 1;
tline = fgetl(fid);
A{i}= tline;
 while ischar(tline)
     i = i+1;
    tline = fgetl(fid);
    A{i} = tline;
  
   end
fclose(fid);

for z = 1:s 
           if matfri + parameters.piper.rotationwidth > 70000
            BlockIndex = results.TransIndex( (z - 1) * parameters.piper.rotationwidth + 1 : end);  
         else
            BlockIndex = results.TransIndex( (z - 1) * parameters.piper.rotationwidth + 1 : parameters.piper.rotationwidth * z);
         end
        K = BlockIndex - n * (z -1) * parameters.piper.rotationwidth;
        Block = [];
for q = 1:pdblength
    p = [parameters.piper.md, pdbls{q}];
    A{21} = p;
    fid = fopen('piper_trypsin_read_deb', 'r+');
    
    for I = 1:24
          fprintf(fid, '%s\n', A{I});
    end
    fclose(fid);
    fprintf('matfri = %d Before------------------------- \n',matfri);
%     matfri=0;
%     EnergiesVec = [];
    fprintf('matfri = %d After------------------------- \n',matfri);
%     for z = 1 : s
%         save('DebugDataMat','-v7.3');
%         fprintf('q =  %d; z = %d, matfri = %d ------------------------- \n',q,z,matfri);
%         [Energies,nx,ny,nz]= testpiperreadDebug(matfri,parameters);
%          if matfri + parameters.piper.rotationwidth > 70000;
%             BlockIndex = results.TransIndex( (z - 1) * parameters.piper.rotationwidth + 1 : end);  
%          else
%             BlockIndex = results.TransIndex( (z - 1) * parameters.piper.rotationwidth + 1 : parameters.piper.rotationwidth * z);
%          end
%         K = BlockIndex - n * (z -1) * parameters.piper.rotationwidth;
%         BlockEnergies = Energies(K);
%         clear Energies;
%         matfri = z * parameters.piper.rotationwidth;
%         EnergiesVec = [EnergiesVec; BlockEnergies];
%      
%     end

%         save('DebugDataMat','-v7.3');
        fprintf('q =  %d; z = %d, matfri = %d ------------------------- \n',q,z,matfri);
        [Energies,nx,ny,nz]= testpiperreadDebug(matfri,parameters);
  
        BlockEnergies = Energies(K);
        clear Energies;
%         matfri = z * parameters.piper.rotationwidth;
%         EnergiesVec = [EnergiesVec; BlockEnergies];
        Block = [Block BlockEnergies];
        save('DebugDataMat','-v7.3');
   
   
end
 matfri = z * parameters.piper.rotationwidth + 36222;
 EnergiesMat = [EnergiesMat; Block];
end
[tf, Index_order]= ismember(results.OrigRotation,results.SortedRotation);
% [results.OrigRotation, Index_order] = sort(results.SortedRotation);
EnergiesMat = EnergiesMat(Index_order,:);

 if SaveEnergy == true
        save(filename, 'EnergiesMat','-v7.3');
 end
    
end
results.EnergiesMat = EnergiesMat;

 allMin = min(results.EnergiesMat, [],'all');
  results.EnergyMin = allMin;
 if results.EnergyMin < 1
    results.ScaledEnergies = EnergiesMat + ( 1 - results.EnergyMin);
    else
    results.ScaledEnergies = EnergiesMat;
 end
    logEnergies = log(results.ScaledEnergies);
results.var = var(logEnergies');
results.mean = mean(logEnergies');

save(filename, 'EnergiesMat','-v7.3');


B(:,12) = results.var;
B(:,13) = results.mean;
results.TransMat = B;
Transfname = append('TransMat','_','dim',num2str(parameters.GRFModel.numofeigenvalues),'_deb');
if parameters.flags.saveTrans == true

save (Transfname, 'B','-v7.3');
end
end