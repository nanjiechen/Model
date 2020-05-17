function [parameters,data] = CollectLocations(methods,parameters);


% Initialize piper
parameters = methods.piperinit(parameters);

% Dummy 
data = [];
%list all the rec_grid files
Savedata = int16(1);
Loaddata = int16(1);
%load file if exists
if exist('xyz.mat', 'file') == 2 & Loaddata == int16(1)
    load('xyz.mat');
else

    [parameters,data] = methods.piper.writedx(methods,parameters);
    %Save data 
    %Flag True for saving data
    if Savedata == int16(1)
        save('xyz.mat', 'data','-v7.3');
    end
    
end

end