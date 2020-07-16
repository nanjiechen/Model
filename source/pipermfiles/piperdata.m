function [parameters,data] = piperdata(methods,parameters);


% Initialize piper
parameters = methods.piperinit(parameters);

% Dummy 
data = [];
%list all the rec_grid files
Savedata = int16(1);
%load file if exists
if exist('dxmatrices.mat', 'file') == 2
    load('dxmatrices.mat');
else

    [parameters,data] = methods.piper.writedx(methods,parameters);
    %Save data 
    %Flag True for saving data
    if Savedata == int16(1)
        save('dxmatrices.mat', 'data','-v7.3');
    end
    
end
end




