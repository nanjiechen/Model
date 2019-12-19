function [parameters,data] = piperdata(methods,parameters);


% Initialize piper
parameters = methods.piperinit(parameters);

% Dummy 
% data = 1;
%list all the rec_grid files
d = '/projectnb/uqproj/StochasticDocking/Code/source/piper/orig_grids';
files = dir(fullfile(d,'rec_grid*'));
ls = {files.name};
Savedata = int16(1);
%load file if exists
if exist('dxmatrices.mat', 'file') == 2
    load('dxmatrices.mat');
else
   % create dx files from list of pdb files
    l = length(ls); 
    %l = 2; %for test
    M = cell(l,1);
    for i = 1:l
        s = ['./piper/orig_grids/' ls{i}];
        parameters.piper.DxFile = s;
        [parameters,localdata] = methods.piper.writedx(methods,parameters);
        M{i} = localdata;
    end
    
    %Save data 
    %Flag True for saving data
    if Savedata == int16(1)
        save('dxmatrices.mat', 'M');
    end
    
end

data.M = cell2mat(M');


