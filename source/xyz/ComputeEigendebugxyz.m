%Construct eigenfunction by using eigenvectors
function [parameters,data] = ComputeEigendebugxyz(parameters,data);

% Compute eigenvectors and eigenvalues of covariance matrix
% m = size(data.protein,1)
% r = data.densityMatrix(:,:,1);

x = data.x; 
y = data.y;
z = data.z;
% mx = size(x,1);
% my = size(y,1);
% mz = size(z,1);
k = parameters.GRFModel.numofeigenvalues;
cx = parameters.GRFModel.cx; 
cy = parameters.GRFModel.cy; 
cz = parameters.GRFModel.cz; 
M = size(cx,2);
% My = size(cy,2);
% Mz = size(cz,2);
Covx = (cx' * cx) / M; %Covariance Matrix
Covy = (cy' * cy) / M;
Covz = (cz' * cz) / M;

% Run SVD to compute eigenvalues and eigenvectors
% [u,s,v] = svds(@(y,tflag) computeevec(y,tflag,cx), [m m], k);
% [u,s,v] = svd[u,s,v] = svds(@(y,tflag) computeevec(y,tflag,cx), [m m], k);s(@(y,tflag) computeevec(y,tflag,cx), [m m], k);
% Mat = [cx;cy;cz];
% C = (Mat' * Mat)/M;
C = Covx + Covy + Covz;
[u,s,v] = svd(C);

s = diag(s);
s = s(1 : k);
u = u(:,1 : k);


% Double check

parameters.GRFModel.covariance = C;
parameters.GRFModel.eigenvalues = s;
parameters.GRFModel.eigenvectors = u;


xefun = parameters.GRFModel.eigenvectors' * cx' ./s /M;
% xefun = cx * u;
% xefun = parameters.GRFModel.eigenvectors' * cx';
parameters.GRFModel.xeigenfunction = xefun;

yefun = parameters.GRFModel.eigenvectors' * cy' ./s /M;
% yefun = cy * u;
% yefun = parameters.GRFModel.eigenvectors' * cy';
parameters.GRFModel.yeigenfunction = yefun;

zefun = parameters.GRFModel.eigenvectors' * cz' ./s /M;
% zefun = cz * u;
% zefun = parameters.GRFModel.eigenvectors' * cz';
parameters.GRFModel.zeigenfunction = zefun;






end