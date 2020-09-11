
% Test comparison between EPFL Sparse Grid Matlab code
% and Kimkle Sparse Grid Code

clear all;
%clc;
%f = @(x,b) prod(1./sqrt(x+b));
%f = @(x,b) prod(exp(-1/2 * x.^2 )) / (sqrt(2*pi)^(size(x,1)));

f = @(x,b) prod(1./ ( (x.^2 + b)));
h = @(x,b) prod(exp(-1/2 * x.^2 )) / (sqrt(2*pi)^(size(x,1)));
g = @(x,b) prod(1./ ( (x.^2 + b))) .* prod(exp(-1/2 * x.^2 )) / (sqrt(2*pi)^(size(x,1)));

b = 3;
N = 2;
I_1d=(2*sqrt(1+b)-2*sqrt(-1+b));
I_ex = I_1d^N;
width = sqrt(3);
width = 6;

% generate the knots and the SM grid. 'nonprob' means we are integrating w.r.t. the pdf rho(x)=1 and not rho(x)=1/prod(b_i - a_i)
%knots=@(n) knots_CC(n,-width,width,'nonprob');
knots=@(n) knots_gaussian(n,0,1);
w = 5;
S = smolyak_grid(N,w,knots,@lev2knots_doubling);
Sr = reduce_sparse_grid(S);

% compute integral
I=f([Sr.knots],b)*[Sr.weights]'  %#ok<NOPTS>




%% Kimkle Sparse Grid
% options = spset('FunctionArgType', 'vector','MinDepth',w,'MaxDepth',w, 'GridType','Chebyshev');
% range = (width*[-1 1]);
% ranges = repmat(range,N,1);
% SG = spvals(@(x) g(x,b), N, ranges, options);
% IKSG = spquad(SG);


%% CC quadrature
w = 15;
knots=@(n) knots_CC(n,-width,width,'nonprob');
S = smolyak_grid(N,w,knots,@lev2knots_doubling);
Sr = reduce_sparse_grid(S);
 
ICC=g([Sr.knots],b)*[Sr.weights]'  %#ok<NOPTS>
 
 
disp('----------');
disp('difference between values');

abs(I-ICC) / abs(ICC)  %#ok<MNEFF,NOPTS>

% compare with exact value
%disp('----------');
%disp('quad error');
%abs(I-I_ex)
