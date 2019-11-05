function [parameters,data] = piperdata(methods,parameters);


% Initialize piper
parameters = methods.piperinit(parameters);

% Dummy 
data = 1;

%load file if exists
if exist('dxdata.mat', 'file') == 2
    load('dxdata.mat');
end

% create dx files from list of pdb files.
[parameters,data] = methods.piper.writedx(methods,parameters);

%Save data 
Savedata = int16(1); %Flag True for saving data
if Savedata == int16(1)
   save('dxdata.mat', 'data');
end

% unload library