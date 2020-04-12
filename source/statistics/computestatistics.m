function results = computestatistics(parameters, methods, results);

% Compute mean and std of the voltages
Sr = parameters.stochasticmodel.Sr;
width = parameters.stochasticmodel.width;

% Collect data
meanData = results.Energes;
VarData  = meanData.^2;

% compute integral
Qmean        = meanData * [Sr.weights]';
Qmeansquare  = VarData  * [Sr.weights]';
    
% Statistics
results.mean = Qmean / ( (2 * width)^parameters.stochasticmodel.N);
    
% Sparse grid mean square
meansquare = Qmeansquare / ( (2 * width)^parameters.stochasticmodel.N);
    
results.var = results.mean - meansquare;