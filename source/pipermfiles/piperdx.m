function [parameters,data] = piperdx(methods,parameters);
% run('paths.m');
% List of pdb files
d = '/projectnb/uqproj/StochasticDocking/Code/data/md_trypsin/';
%files = dir(['/projectnb/uqproj/StochasticDocking/Code/data/md_trypsin/','*.pdb']);
files = dir([d,'*.pdb']);
%files =dir([parameters.piper.md,'/*.pdb']);

list = {files.name};
%l = length(list);
l = 1; %For test

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
data = [];
for j = 1:l
    p = append(parameters.piper.md,list{j});
    A{19} = p;

   
    fid = fopen(parameters.piper.piperfile, 'r+');
    %for I = 1:length(A)
    for I = 1:23
        fprintf(fid, '%s\n', char(A{I}));
    end
    fclose(fid);
    run('testpipergridonly.m');
    d = dx2mat(parameters.piper.DxFile);
    dmatrix = d.densityMatrix;
    dvector = dmatrix(:);
    data = [data dvector];
end

%data = dx2mat(parameters.piper.DxFile);


end



% for each run piper

% read the dx file