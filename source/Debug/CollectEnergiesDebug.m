
%Obtain Energies
function [parameters,results] = CollectEnergiesDebug(matfri,parameters)
parameters.piper.rotationwidth = 3000;
matfri = 0;
% l = parameters.stochasticmodel.numofgrids;
n = parameters.piper.nx * parameters.piper.ny * parameters.piper.nz;
% inputd = '/projectnb/uqproj/StochasticDocking/Code/source/piper/input_grids';
% 
EnergiesMat = [];
A = [];
s = ceil(69999/parameters.piper.rotationwidth);
% % List of pdb files
 data = [];
% d = '/projectnb/uqproj/StochasticDocking/Code/data/md_trypsin/';
% pdbfiles = dir([d,'*out*.pdb']);
% pdbls = {pdbfiles.name};
% pdblength = length(pdbfiles);

SaveEnergy = parameters.flags.saveenergies;
LoadEnergy = parameters.flags.loadenergies;

%filename = 'energiesMat.mat';
filename = append('energiesMat','_','dim',num2str(parameters.GRFModel.numofeigenvalues),'_deb','.mat');
if exist(filename, 'file') == 2 & LoadEnergy == true
    load(filename);
else

for z = 1:s  
    matfri = (z -1) * parameters.piper.rotationwidth;

    fprintf('matfri = %d ------------------------- \n',matfri);

  [Energies,nx,ny,nz]= testpiperreadDebug(matfri,parameters);
  matlri = matfri + parameters.piper.rotationwidth - 1;
  b = [];
  if matlri <= 69999
      w = parameters.piper.rotationwidth;
  else
      w = 70000 - matfri;
  end
  for r = 1:w
      energyvec = Energies(n*(r-1)+1 : n*r);
      b(1) = r + (z-1)*parameters.piper.rotationwidth;
      [e,f] = min(energyvec);
      b(2) = f;
      b(3) = e;
      A = [A;b];
  end
        
        clear Energies;
        save('DebugDataMat','-v7.3');
   
   
end

end

end
