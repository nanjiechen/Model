function [parameters,results] = EnergiesBlocks(parameters,results)
n = ceil(69999/parameters.piper.rotationwidth);
matfri = 0;
%filename = append('energiesMat','_','dim',num2str(parameters.GRFModel.numofeigenvalues),'.mat');
%mf = matfile(filename, 'Writable', true);
%[nrows,ncols] = size(mf,'EnergiesMat');
%u = nrows/parameters.piper.nx/parameters.piper.ny/parameters.piper.nz/parameters.piper.rotationwidth;
%matfri = parameters.piper.rotationwidth * u;
%k = 1000;
Emin = [];
EVar = [];
Minsindex = [];
nx = parameters.piper.nx;
ny = parameters.piper.ny;
nz = parameters.piper.nz;
Minsrotation = [];
%TransMat = [];
 
[parameters,results] = CollectEnergies(matfri,parameters,results);

 eg = results.EnergiesMat;
 allMin = min(results.EnergiesMat, [],'all');
  results.EnergyMin = allMin;
 if results.EnergyMin < 1
    results.ScaledEnergies = eg + ( 1 - results.EnergyMin);
    else
    results.ScaledEnergies = eg;
 end
   
 
   
    %((i - 1)* parameters.piper.rotationwidth + 1 : i * parameters.piper.rotationwidth, :)
   
    [parameters,results] = ComputeStats(parameters,results);
   
   


%Ordered back to the unsorted original order.
[tf, Index_order]= ismember(results.OrigRotation,results.SortedRotation);
% [results.OrigRotation, Index_order] = sort(results.SortedRotation);
ReorderedMat = results.TransMat(Index_order,:);
results.ReorderedMat = ReorderedMat;
Reorderedfname = append('ReorderedMat','_','dim',num2str(parameters.GRFModel.numofeigenvalues));
if parameters.flags.saveReordered == true
    save(Reorderedfname,'ReorderedMat','-v7.3');
end

end
