% Build a map between sparse grids and K-L approximation
function [parameters] = RecMap(methods,parameters)
    parameters = methods.sparseinit(parameters);
    knots = parameters.stochasticmodel.Sr.knots;
    evals = parameters.GRFModel.eigenvalues;
    evals = sqrt(evals);
    efunc = parameters.GRFModel.eigenfunction;
    d = diag(evals);
    prod = efunc'*d;
    l = size(parameters.stochasticmodel.Sr.m);
    m = l(1); %Number of grids
    parameters.stochasticmodel.numofgrids = m;
    MeanMat = repmat(transpose(parameters.GRFModel.meanvector),1,m);
    parameters.stochasticmodel.RecMap = prod*knots + MeanMat ;
end