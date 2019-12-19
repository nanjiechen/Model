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

%Create Sparse grids
methods.sparseinit = @sparseinitialization;

%Karhunen-Loeve expansion for receptor map.
methods.RecModel = @RecMap;

% Define parameters for GRF model
parameters.GRFModel.numofeigenvalues = 2;

% Load stochastic realizations
[parameters, data] = methods.loadrealizations(methods,parameters);

% Create stochastic protein model
[parameters,data] = methods.ProteinStochasticModel(methods,parameters,data);
%Get the approximate receptor map by using Karhunen-Loeve expansion
[parameters] = RecMap(methods,parameters);


