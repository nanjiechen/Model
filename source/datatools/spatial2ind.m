function [parameters,results] = spatial2ind(parameters)

% Load Piper Results Matrix
A = load('../toolbox/dimitritools/ft.000.00');
M = 70000;
% fft_tind_to_gind.pl 0.79   0.46   0.79 84 96 84 1.0 13.792500 -1.544000 11.785500
if parameters.loadTransindex ==  true & exist('spatial2ind.mat', 'file')
      load('spatial2ind.mat');
else
for i = 1:M
    command = ['../toolbox/dimitritools/fft_tind_to_gind.pl ',num2str(A(i,2:4)), ' 96 90 80 1.0 20.289500 3.523500 -5.085500'];

    %command = '../toolbox/dimitritools/fft_tind_to_gind.pl 0.79   0.46   0.79 84 96 84 1.0 13.792500 -1.544000 11.785500';
    [status,cmdout] = system(command);

    a = splitlines(cmdout);

    transindex = str2num(a{3}) + 1;

    A(i,11) = transindex;
    
end
end
save('spatial2ind.mat','A');
transition = A(:,11);
% OrigInd = transition + A(:,1)* parameters.piper.nx * parameters.piper.ny * parameters.piper.nz;
% results.OrigInd = OrigInd;
[B,I] = sort(A(:,1));
results.OrigRotation = A(:,1);
results.SortedRotation = B;
%SI is the sorted transition
SI = transition(I);
TransInd = SI + B * parameters.piper.nx * parameters.piper.ny * parameters.piper.nz;
results.spatial2ind = A;
results.Index = I;
results.TransIndex = TransInd;
% J = results.TransIndex <= parameters.piper.nx * parameters.piper.ny * parameters.piper.nz * parameters.piper.rotationwidth;
% results.NewIndex = results.TransIndex(J);
%results.NewIndex = TransInd;
AscendMat = A(I,:);
results.SortedMat = AscendMat;

end