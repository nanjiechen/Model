function [parameters,data] = ComputeEigen(parameters,data);

%m = size(data.protein,1)
%r = data.densityMatrix(:,:,1);
%r = data;
r = data.M; % dummy. Needs to be corrected.

m = size(r,1);
k = parameters.GRFModel.numofeigenvalues;
cx = parameters.GRFModel.cx; %Covariance Matrix

% Run SVD to compute eigenvalues and eigenvectors
[u,s,v] = svds(@(y,tflag) computeevec(y,tflag,cx), [m m], k);


parameters.GRFModel.eigenvalues = diag(s);
parameters.GRFModel.eigenvectors = u;


function y = computeevec(x,tflag,cx); 
    y = cx * x;
    y = cx' * y; 
end

end