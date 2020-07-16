function parameters = sparseinitialization(parameters)
    parameters.stochasticmodel.N = parameters.GRFModel.numofeigenvalues;
    parameters.stochasticmodel.width = 3;
    parameters.stochasticmodel.level = 3;
    methods.sparsegridinit = @sginit;
    parameters = sginit(parameters);
end