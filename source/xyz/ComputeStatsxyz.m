function [parameters,results] = ComputeStatsxyz(parameters,results)
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
varlogData = meanlogData.^2;
results.SampleMean = mean(meanlogData');
results.SampleVariance = var(meanlogData');
% compute integral


% Qmeanlog = [meanlogData .* parameters.GRFModel.density(knots)] * [Sr.weights]';
% Qmeanlogsquare  = [varlogData .* parameters.GRFModel.density(knots)] * [Sr.weights]';
Qmeanlog = [meanlogData] * [Sr.weights]';
Qmeanlogsquare  = [varlogData] * [Sr.weights]';
% Statistics
results.meanlog = Qmeanlog;

% Sparse grid mean square
results.meansquare = Qmeanlogsquare;
    
results.varlog = results.meansquare - results.meanlog.^2;

A = results.SortedMat;

A(:,12) = results.varlog;
A(:,13) = results.meanlog;

results.TransMat = A;
Transfname = append('../data/PiperData/','TransMat','_','dim',num2str(parameters.GRFModel.numofeigenvalues),'xyz');
if parameters.flags.saveTrans == true

save (Transfname, 'A','-v7.3');
end
save '../data/PiperData/PiperData';
end