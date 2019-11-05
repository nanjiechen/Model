% Methods:

%1] Method to load data from realizations MD dynamic simulations

%2] Code to contruct the Gaussian Random Field model of the nprotein

%3] Construct optimization of sochastic docking

%4] Perform stochastic docking

%5] Output results

% Define methods


% Piper
methods.piperinit = @piperinitialization;

% Data types
% methods.loadrealizations = @fakeproteinrealizations;

% Piper methods
methods.loadrealizations = @piperdata;
methods.piper.writedx = @piperdx;

% Stochastic model methods
methods.ProteinStochasticModel = @GRFmodel;

% Eigenmodels
methods.EigenDecomposition = @ComputeEigen;

% Define parameters for GRF model
parameters.GRFModel.numofeigenvalues = 3;

% Load stochastic realizations
[parameters, data] = methods.loadrealizations(methods,parameters);

% Create stochastic protein model
[parameters, data] = methods.ProteinStochasticModel(methods,parameters,data);