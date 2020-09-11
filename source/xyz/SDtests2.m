% Methods:

%1] Method to load data from realizations MD dynamic simulations

%2] Code to contruct the Gaussian Random Field model of the nprotein

%3] Construct optimization of sochastic docking

%4] Perform stochastic docking

%5] Output results

% Define methods



% Data types
% methods.loadrealizations = @fakeproteinrealizations;
% Piper1
methods.piperinit = @piperinitializationxyz;
% Piper methods
methods.loadrealizations = @CollectLocations;
methods.piper.writedx = @piperxyz;

% Stochastic model methods
methods.ProteinStochasticModel = @GRFmodelxyz;

% Eigenmodels
methods.EigenDecomposition = @ComputeEigenxyz;

% Create Sparse grids
methods.sparseinit = @sparseinitialization;

% Karhunen-Loeve expansion for receptor map.
methods.RecModel = @RecMapxyz;
% methods.Spatial2ind = @spatial2ind;
 methods.Spatial2ind = @SortEnergies;
% Run Piper again, get engergies.
% methods.Energies = @EnergiesBlocks;
methods.writepdb = @writepdb;
methods.goodmultirec = @multirec;

methods.CollectEng = @CollectEnergiesxyz;
methods.ComputeStats = @ComputeStatsxyz;

% Define parameters for GRF model
parameters.GRFModel.numofeigenvalues = 2;
parameters.GRFModel.KLRescale = 0.001;

%Define flag parameters
parameters.flags.saveenergies = true;
parameters.flags.loadenergies = true;
parameters.flags.saveTrans = true;
parameters.loadTransindex = true;
parameters.flags.saveReordered = true;
% memory parameters
parameters.piper.memory =   1.8e+10;
%parameters.piper.memory =   4.8e+10;




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
% [parameters,approc] = debug(parameters);
fprintf('Load Energies and Compute statistics -------------------- \n'); 
fprintf('\n');
% Obtain Energies
[parameters] = methods.writepdb(parameters);
[parameters] = methods.goodmultirec(parameters);
% Initialize piper
parameters = methods.piperinit(parameters);
[parameters,results] = methods.Spatial2ind(parameters);
% [parameters,results] = methods.CollectEng(parameters,results);
[parameters,results] = methods.CollectEng(parameters,results);
% Compute Statistics
fprintf('Compute statistics -------------------------------------- \n');
fprintf('\n');
[parameters,results] = methods.ComputeStats(parameters,results);

% Unload piper library
unloadlibrary('libMatPiper')