 % Gaussian Random Field model from data

function [parameters,data] = GRFmodelxyz(methods,parameters,data);

%r = data.protein;
%r = parameters.densityMatrix;
%r = data.densityMatrix(:,:,1);

x = data.x;
y = data.y;
z = data.z;
mx = size(x,2);
my = size(y,2);
mz = size(z,2);
mat = data.mat;
meanx = mean(x');
xmeanmatrix = repmat(transpose(meanx),1,mx);
parameters.GRFModel.xmeanmatrix = xmeanmatrix;
parameters.GRFModel.xmeanvector = meanx;
cx = x - xmeanmatrix;

meany = mean(y');
ymeanmatrix = repmat(transpose(meany),1,my);
parameters.GRFModel.ymeanmatrix = ymeanmatrix;
parameters.GRFModel.ymeanvector = meany;
cy = y - ymeanmatrix;

meanz = mean(z');
zmeanmatrix = repmat(transpose(meanz),1,mz);
parameters.GRFModel.zmeanmatrix = zmeanmatrix;
parameters.GRFModel.zmeanvector = meanz;
cz = z - zmeanmatrix;

% Density model
cte = (1 / (2 * pi)^(parameters.GRFModel.numofeigenvalues/2));
parameters.GRFModel.density = @(knots) prod(exp(-1/2 * (knots).^2 )) * cte;

parameters.GRFModel.cx = cx;
parameters.GRFModel.cy = cy;
parameters.GRFModel.cz = cz;
[parameters,data] = methods.EigenDecomposition(parameters,data);
end