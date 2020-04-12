function [parameters,results] = EnergiesBlocks(parameters,results)
%n = ceil(69999/parameters.piper.rotationwidth);
matfri = 0;
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

 eg = results.EnergiesMat(results.NewIndex);
 allMin = min(results.EnergiesMat(results.NewIndex), [],'all');
  results.EnergyMin = allMin;
 if results.EnergyMin < 1
    results.ScaledEnergies = eg + ( 1 - results.EnergyMin);
    else
    results.ScaledEnergies = eg;
 end
   
 
   
    %((i - 1)* parameters.piper.rotationwidth + 1 : i * parameters.piper.rotationwidth, :)
   
    [parameters,results] = ComputeStats(parameters,results);
   
   


%Ordered back to the unsorted original order.
[results.UnsortedIndex, Index_order] = sort(results.NewIndex);
ReorderedMat = results.TransMat(Index_order,:);
results.ReorderedMat = ReorderedMat;
Reorderedfname = append('ReorderedMat','_','dim',num2str(parameters.GRFModel.numofeigenvalues));
if parameters.flags.saveReordered == true
    save(Reorderedfname,'ReorderedMat','-v7.3');
end

end
