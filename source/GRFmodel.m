 % Gaussian Random Field model from data

function [parameters,data] = GRFmodel(methods,parameters,data);

%r = data.protein;
%r = parameters.densityMatrix;
%r = data.densityMatrix(:,:,1);

r = data;

m = size(r,2);

meanvector = mean(r');
meanmatrix = repmat(transpose(meanvector),1,m);
parameters.GRFModel.meanmatrix = meanmatrix;
parameters.GRFModel.meanvector = meanvector;
cx = r - meanmatrix;

% Density model
cte = (1 / (2 * pi)^(parameters.GRFModel.numofeigenvalues/2));
parameters.GRFModel.density = @(knots) prod(exp(-1/2 * (knots).^2 )) * cte;

parameters.GRFModel.cx = cx;

[parameters,data] = methods.EigenDecomposition(parameters,data);
end