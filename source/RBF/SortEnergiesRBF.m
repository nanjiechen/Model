%Obtain Energies
function [parameters,results] = SortEnergiesRBF(parameters)
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

if exist('SortedMatrix.mat', 'file') == 2 
    load('SortedMatrix.mat');
else

for z = 1:s  
    matfri = (z -1) * parameters.piper.rotationwidth;

    fprintf('matfri = %d ------------------------- \n',matfri);

  [Energies,nx,ny,nz]= testpiperreadSorting(matfri,parameters);
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
        save('Data\Sorthing.mat','-v7.3');
   
   
end

end

[S,ind] = sort(A(:,3));
A = A(ind,:);
save('../data/PiperData/SortedMatrix.mat','A');
results.OrderedMatrix = A;
%  v = [802,3118,4832,7933,8657,8907,16771,8502,20448,24890]'; 
 v = [802]';
D = A(v,:);
SortedMat = D;
results.SortedMat = D;
%% 
 z = parameters.piper.nx * parameters.piper.ny * parameters.piper.nz;
% load('SortedMat.mat');
% transition = SortedMat(:,2);
% [B,I] = sort(SortedMat(:,1));
% results.OrigRotation = SortedMat(:,1);
% results.SortedRotation = B;
% SI = transition(I);
% TransInd = SI + (B-1) * z;
% results.spatial2ind = SortedMat;
% results.Index = I;
% results.TransIndex = TransInd;
% AscendMat = SortedMat(I,:);
% results.SortedMat = AscendMat;
results.TransIndex = SortedMat(:,2) + z * (SortedMat(:,1) -1);
%% 
end

