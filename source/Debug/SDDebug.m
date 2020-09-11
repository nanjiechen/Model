% Methods:

%1] Method to load data from realizations MD dynamic simulations

%2] Code to contruct the Gaussian Random Field model of the nprotein

%3] Construct optimization of sochastic docking

%4] Perform stochastic docking

%5] Output results

% Define methods

% Piper1
methods.piperinit = @piperinitializationDebug;

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
methods.Spatial2ind = @spatial2ind;
% Run Piper again, get engergies.

methods.CollectEng = @CollectEnergiesDebug;


% Define parameters for GRF model
parameters.GRFModel.numofeigenvalues = 2;
parameters.flags.saveenergies = true;
parameters.flags.loadenergies = false;
parameters.flags.saveTrans = true;
parameters.loadTransindex = true;
parameters.flags.saveReordered = true;
matfri = 0;
% memory parameters
parameters.piper.memory = 2.0e+11;

% Diary 
diary my_simulation.out
% Load stochastic realizations using piper
fprintf('Load realizations --------------------------------------- \n');
fprintf('\n');
parameters = methods.piperinit(parameters);
% Load stochastic realizations using piper
% fprintf('Load realizations --------------------------------------- \n');
% fprintf('\n');
% [parameters, data] = methods.loadrealizations(methods,parameters);

% fprintf('Create stochastic protein model ------------------------- \n');
% fprintf('\n');
% 
% % Create stochastic protein model
% [parameters,data] = methods.ProteinStochasticModel(methods,parameters,data);
% 
% fprintf('Stochastic model ---------------------------------------- \n');
% fprintf('\n');
% % Get the approximate receptor map by using Karhunen-Loeve expansion
% [parameters] = methods.RecModel(methods,parameters);
% 
% fprintf('Load Energies and Compute statistics -------------------- \n'); 
% fprintf('\n');
% % Obtain Energies
% [parameters,results] = methods.Spatial2ind(parameters);
[parameters,results] = methods.CollectEng(matfri,parameters);

% Compute Statistics
fprintf('Compute statistics -------------------------------------- \n');
fprintf('\n');
%[parameters,results] = methods.ComputeStats(parameters,results);

% Unload piper library
unloadlibrary('libMatPiper')

diary off;