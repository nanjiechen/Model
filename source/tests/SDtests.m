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

% Create Sparse grids
methods.sparseinit = @sparseinitialization;

% Karhunen-Loeve expansion for receptor map.
methods.RecModel = @RecMap;

% Run Piper again, get engergies.
methods.Energies = @CollectEnergies;
methods.ComputeStats = @ComputeStats;

% Define parameters for GRF model
parameters.GRFModel.numofeigenvalues = 2;
parameters.flags.saveenergies = true;
parameters.flags.loadenergies = true;



% Load stochastic realizations using piper
fprintf('Load realizations --------------------------------------- \n');
fprintf('\n');
[parameters, data] = methods.loadrealizations(methods,parameters);

fprintf('Create stochastic protein model ------------------------- \n');
fprintf('\n');

% Create stochastic protein model
[parameters,data] = methods.ProteinStochasticModel(methods,parameters,data);

fprintf('Stochastic model ---------------------------------------- \n');
fprintf('\n');
% Get the approximate receptor map by using Karhunen-Loeve expansion
[parameters] = methods.RecModel(methods,parameters);


fprintf('Load Energies ------------------------------------------- \n');
fprintf('\n');
% Obtain Engernies
[parameters,results] = methods.Energies(parameters);

% Compute Statistics
fprintf('Compute statistics -------------------------------------- \n');
fprintf('\n');
[parameters,results] = methods.ComputeStats(parameters,results);

% Unload piper library
unloadlibrary('libMatPiper')
