%Construct eigenfunction by using eigenvectors
function [parameters,data] = ComputeEigen(parameters,data);

% Compute eigenvectors and eigenvalues of covariance matrix
% m = size(data.protein,1)
% r = data.densityMatrix(:,:,1);

r = data; 

m = size(r,1);
k = parameters.GRFModel.numofeigenvalues;
cx = parameters.GRFModel.cx; 

M = size(cx,2);

C = (cx' * cx) / M; %Covariance Matrix

% Run SVD to compute eigenvalues and eigenvectors
% [u,s,v] = svds(@(y,tflag) computeevec(y,tflag,cx), [m m], k);
% [u,s,v] = svd[u,s,v] = svds(@(y,tflag) computeevec(y,tflag,cx), [m m], k);s(@(y,tflag) computeevec(y,tflag,cx), [m m], k);

[u,s,v] = svd(C);

s = diag(s);
s = s(1 : k);
u = u(:,1 : k);

% Double check

parameters.GRFModel.eigenvalues = s;
parameters.GRFModel.eigenvectors = u;

efun = parameters.GRFModel.eigenvectors' * cx';
parameters.GRFModel.eigenfunction = efun;
% function y = computeevec(x,tflag,cx); 
%    y = cx * x;
%    y = 1/100 * cx' * y; 
%end

end