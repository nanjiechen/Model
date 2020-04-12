%Obtain Energies
function [parameters,results] = CollectEnergies(parameters)
SaveEngergies = int16(1);
d = '/projectnb/uqproj/StochasticDocking/Code/data/md_trypsin/';
files = dir([d,'*.pdb']);
list = {files.name};
m = parameters.stochasticmodel.Sr.m;
l = length(m);
%l = 2;
EnergiesMat = [];

if exist('Energies.mat', 'file') == 2
    load('Energies.mat');
else
    fid = fopen(parameters.piper.piperfile,'r');
    i = 1;
    tline = fgetl(fid);
    A{i}= tline;
    while ischar(tline)
     i = i+1;
     tline = fgetl(fid);
     A{i} = tline;
  
   end
fclose(fid);

    for i = 1:l
    p = [parameters.piper.md, list{i}];
    A{19} = p;
       fid = fopen(parameters.piper.piperfile, 'r+');
       for I = 1:22
      
        fprintf(fid, '%s\n', A{I});
    end
    fclose(fid);
    run('testpiper.m');
    EnergiesMat = [EnergiesMat Energies];
    end
    
    if SaveEngergies== int16(1)
        
       save('Energies.mat', 'EnergiesMat','-v7.3');
     end

end
results.energies = EnergiesMat;

end