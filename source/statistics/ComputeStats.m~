function [parameters,results] = ComputeStats(parameters,results)
% Compute mean and std of the voltages
Sr = parameters.stochasticmodel.Sr;
width = parameters.stochasticmodel.width;
knots = parameters.stochasticmodel.Sr.knots;
n = parameters.stochasticmodel.N;

% x1 = knots(1,:);
% x2 = knots(2,:);

% Collect data
meanData = results.ScaledEnergies; 
VarData  = meanData.^2;
meanlogData = log(meanData);
% compute integral

Qmean        = [meanData .* [Sr.weights]] * parameters.GRFModel.density(knots)';
Qmeansquare  = [VarData .* [Sr.weights]] * parameters.GRFModel.density(knots)';
Qmeanlog = [meanlogData .* [Sr.weights]] * parameters.GRFModel.density(knots)';

% Statistics
results.mean = Qmean / ( (2 * width)^parameters.GRFModel.numofeigenvalues);
results.meanlog = Qmeanlog /( (2 * width)^parameters.GRFModel.numofeigenvalues);
% Sparse grid mean square
meansquare = Qmeansquare / ( (2 * width)^parameters.GRFModel.numofeigenvalues);
    
results.var = meansquare - results.mean.^2;

A = results.SortedRotation;

A(:,12) = results.var;
A(:,13) = results.mean;

results.TransMat = A;
Transfname = append('TransMat','_','dim',num2str(parameters.GRFModel.numofeigenvalues));
if parameters.flags.saveTrans == true

save (Transfname, 'A','-v7.3');
end
end