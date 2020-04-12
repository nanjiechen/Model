function results = computestatisticsperlevel(parameters, methods, results);

% Compute mean and std of the voltages
Sr = parameters.stochasticmodel.Sr;
width = parameters.stochasticmodel.width;
level = parameters.stochasticmodel.level;
Xorg = [Sr.knots]';
knots=@(n) knots_CC(n,-width,width,'nonprob');

% Collect data
meanData = results.power.voltages;
VarData  = meanData.^2;

for m = 1 : level
        
    S  =  smolyak_grid(parameters.stochasticmodel.N,m,knots,@lev2knots_doubling);
    Sr = reduce_sparse_grid(S);
    n  = size(Sr.knots,2);
    numofknots(m) = n;
    X = [Sr.knots]';
    [dummy,loc] = ismembertol(X,Xorg,'ByRows',1e-16);
    
    
    % compute integral
    Qmean        = meanData(:,loc) * [Sr.weights]';
    Qmeansquare  = VarData(:,loc) * [Sr.weights]';
    
    % Statistics
    SGmean(:,m) = Qmean / ( (2 * width)^parameters.stochasticmodel.N);
    
    % Sparse grid mean square
    EQmeanSquare(:,m) = Qmeansquare / ( (2 * width)^parameters.stochasticmodel.N);
    
    SGvar(:,m) = EQmeanSquare(:,m) - SGmean(:,m).^2;
end

results.SGmean = SGmean;
results.EQmeanSquare = EQmeanSquare;
results.SGvar = SGvar;
results.numofknots = numofknots;
