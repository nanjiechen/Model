
function parameters = RBFinitialization(parameters)
%% 1. Obtain scatter points to interpolate

    % 0. Obtain the interpolating points, which should be randomly selected points in a 3-std-circle, stored in gridnodes
    % n  Number of random sample points
    n = 10;
    parameters.stochasticmodel.NumOfRBFNodes = n;
        A= 6*rand(n,2)-ones(n,1)*[3,3];
        dist=sqrt(sum(A.^2,2));
        c=(dist<3);
        gridnodes=A(c,:); % (column vector)
        parameters.stochasticmodel.RBF.nodes = gridnodes';
end