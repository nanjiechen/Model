function parameters = sginit(parameters);

% Build sparse grid from 
width = parameters.stochasticmodel.width;
level = parameters.stochasticmodel.level;
N = parameters.stochasticmodel.N;

% generate the knots of the Smolyak (SM) grid
knots=@(n) knots_CC(n,-width,width,'nonprob');
S = smolyak_grid(N,level,knots,@lev2knots_doubling);
parameters.stochasticmodel.Sr = reduce_sparse_grid(S);


%% Maximum rotation width search
parameters.piper.rotationwidth = floor(parameters.piper.memory / ...
    (parameters.piper.memperrot * size(parameters.stochasticmodel.Sr.knots,2)));

end
