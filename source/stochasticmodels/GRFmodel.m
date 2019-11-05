% Gaussian Random Field model from data

function [parameters,data] = GRFmodel(methods,parameters,data);

%r = data.protein;
%r = parameters.densityMatrix;
%r = data.densityMatrix(:,:,1);
r = data;
m = size(r,1);

meanvector = sum(r)/m;
meanmatrix = repmat(transpose(meanvector),1,m);
cx = transpose(r) - meanmatrix;

parameters.GRFModel.cx = cx;

[parameters,data] = methods.EigenDecomposition(parameters,data);
end