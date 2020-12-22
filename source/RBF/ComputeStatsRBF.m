function [parameters,results] = ComputeStatsRBF(parameters,results)
% Compute mean and std of the voltages
Sr = parameters.stochasticmodel.Sr;
width = parameters.stochasticmodel.width;
knots = parameters.stochasticmodel.Sr.knots;
n = parameters.stochasticmodel.N;

% x1 = knots(1,:);
% x2 = knots(2,:);

% Collect data
meanData = results.ScaledEnergies; 
VarData  = meanData.^2;
meanlogData = log(meanData);
varlogData = meanlogData.^2;
results.SampleMean = mean(meanlogData');
results.SampleVariance = var(meanlogData');
% compute integral


% Qmeanlog = [meanlogData .* parameters.GRFModel.density(knots)] * [Sr.weights]';
% Qmeanlogsquare  = [varlogData .* parameters.GRFModel.density(knots)] * [Sr.weights]';
Qmeanlog = [meanlogData] * [Sr.weights]';
Qmeanlogsquare  = [varlogData] * [Sr.weights]';
% Statistics
results.meanlog = Qmeanlog;

% Sparse grid mean square
results.meansquare = Qmeanlogsquare;
    
results.varlog = results.meansquare - results.meanlog.^2;

A = results.SortedMat;

A(:,12) = results.varlog;
A(:,13) = results.meanlog;

results.TransMat = A;
Transfname = append('../data/PiperData/','TransMat','_','dim',num2str(parameters.GRFModel.numofeigenvalues),'RBF');
if parameters.flags.saveTrans == true

save (Transfname, 'A','-v7.3');
end
results.GQStats = [results.meanlog results.varlog];
%% 
s = size(results.SortedMat,1);
SimpsonStats = zeros(s,4);
localmin_vector = zeros(1,s);


for i = 1:s
for j = 1:2    
R = results.R;
R.RBFcoeff = results.RBFCoeffVec(:,i);
% gridnodes = parameters.stochasticmodel.RBF.nodes';
% lowerbounds=[-3, -3];
% upperbounds=[3, 3];
% POLYopts.varnames = {'X','Y'};
% POLYopts.basistype = 'legendre';
% POLYopts.degrees=[4,4]; % degree of polynomials in each dimension
% POLYopts.indexsettype='totaldegree';
% POLYopts.totaldegreevalue = 5;
% RBFopts.shapeparam=[5, 1, 1];
% RBFopts.kerneltype = 5;
    opts.quadtype = 'section';
    opts.numofsections_1d = 100;
    opts.power = 1;
    opts.quadrule = 'simpson'; % rectangle  trapezoidal simpson
    opts.takelog='no';
    [Vol,localmin]= R.Quadrature(opts);
    localmin_vector(i,1)=localmin; 


   opts.globalmin=min(localmin_vector);
   opts.takelog='yes';
   [Vol]= R.Quadrature(opts);
   SimpsonStats(i,j) = Vol;
end

end
results.SimpsonStats = SimpsonStats;
%% 

save '../data/PiperData/PiperResults';
end