function [parameters,data] = piperxyz(methods,parameters);

% List of pdb files
data = [];
d = '/projectnb/uqproj/StochasticDocking/Code/data/md_trypsin/';
%files = dir(['/projectnb/uqproj/StochasticDocking/Code/data/md_trypsin/','*.pdb']);
f = dir([d,'*out*.pdb']);
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
x = [];
y = [];
z = [];
j =1;
for j = 1:l
   
  PDBStruct = pdbread(list{j});
  
  nx = {PDBStruct.Model.Atom.X};
  ny = {PDBStruct.Model.Atom.Y};
  nz = {PDBStruct.Model.Atom.Z};
   save('xyzMat','-v7.3');
  isEmx = cellfun(@isempty, nx);
  nx = cell2mat(nx);
  ny = cell2mat(ny);
  nz = cell2mat(nz);
  nx = nx';
  ny = ny';
  nz = nz';
  nx = single(nx);
  ny = single(ny);
  nz = single(nz);
    
    x = [x nx];
    y = [y ny];
    z = [z nz];
 
end

data.x = x;
data.y = y;
data.z = z;

mat = [x y z];
data.mat = mat;

end