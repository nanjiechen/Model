% Build a map between sparse grids and K-L approximation
function [parameters] = RecMapRBF(methods,parameters)
%     parameters = methods.sparseinit(parameters);
%     knots = parameters.stochasticmodel.Sr.knots;
knots =  parameters.stochasticmodel.RBF.nodes;
    evals = parameters.GRFModel.eigenvalues;
%     evals = sqrt(evals);
    evals = sqrt(evals*parameters.GRFModel.KLRescale);
    xefunc = parameters.GRFModel.xeigenfunction;
    yefunc = parameters.GRFModel.yeigenfunction;
    zefunc = parameters.GRFModel.zeigenfunction;
    d = diag(evals);
    
    xprod = xefunc'*d;
    yprod = yefunc'*d;
    zprod = zefunc'*d;
    
%     l = size(parameters.stochasticmodel.Sr.m);
%     m = l(1);
    m = size(parameters.stochasticmodel.RBF.nodes,2); %Number of grids
    parameters.stochasticmodel.numofgrids = m;
    xMeanMat = repmat(transpose(parameters.GRFModel.xmeanvector),1,m);
    yMeanMat = repmat(transpose(parameters.GRFModel.ymeanvector),1,m);
    zMeanMat = repmat(transpose(parameters.GRFModel.zmeanvector),1,m);
  
    parameters.stochasticmodel.xRecMap = xprod*knots + xMeanMat ;
    parameters.stochasticmodel.yRecMap = yprod*knots + yMeanMat ;
    parameters.stochasticmodel.zRecMap = zprod*knots + zMeanMat ;
end