function [parameters,data] = piperdx(methods,parameters);

% List of pdb files
data = [];
d = '/projectnb/uqproj/StochasticDocking/Code/data/md_trypsin/';
%files = dir(['/projectnb/uqproj/StochasticDocking/Code/data/md_trypsin/','*.pdb']);
f = dir([d,'*.pdb']);
%files =dir([parameters.piper.md,'/*.pdb']);

list = {f.name};
l = length(list);
%l = 2; %For test
d = '/projectnb/uqproj/StochasticDocking/Code/source/piper/orig_grids';
files = dir(fullfile(d,'rec_grid*'));
%List of rec maps
ls = {files.name};
n = length(ls); %Number of rec_grid maps
parameters.piper.NumofRecMap = n;

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


for j = 1:l
   % p = append(parameters.piper.md,list{j});
   p = [parameters.piper.md, list{j}];
    A{19} = p;

   
    fid = fopen(parameters.piper.piperfile, 'r+');
    
    for I = 1:22
        fprintf(fid, '%s\n', A{I});
    end
    fclose(fid);
   run('piperdxwrite.m');
   RecData = [];
   for i = 1:n
        s = ['./piper/orig_grids/' ls{i}];
        parameters.piper.DxFile = s;
        d = dx2mat(parameters.piper.DxFile);
    dmatrix = d.densityMatrix;
    dvector = dmatrix(:); 
    RecData = [RecData;dvector];
    
    end
    data = [data RecData];
    data = single(data);
end




end
