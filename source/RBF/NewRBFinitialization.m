function parameters = NewRBFinitialization(parameters)
load('gridknots.mat');
parameters.stochasticmodel.RBF.nodes = gridnodes';

end
