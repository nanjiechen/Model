function [parameters,data] = CollectEnergiesxyz(parameters);

% List of pdb files
s = ceil(69999/parameters.piper.rotationwidth);
data = [];
matfri = 0;
d = '/projectnb/uqproj/StochasticDocking/Code/data/md_trypsin/';
%files = dir(['/projectnb/uqproj/StochasticDocking/Code/data/md_trypsin/','*.pdb']);
f = dir([d,'*approximation*.pdb']);
%files =dir([parameters.piper.md,'/*.pdb']);

list = {f.name};

n = parameters.stochasticmodel.numofgrids;
md = '../../data/md_trypsin/';

fid = fopen(parameters.piper.piperfileread,'r');
i = 1;
tline = fgetl(fid);
A{i}= tline;
 while ischar(tline)
     i = i+1;
    tline = fgetl(fid);
    A{i} = tline;
  
   end
fclose(fid);
EnergiesMat = [];
 for k = 1:s
      EnergiesVec = [];
%     matfri = 0;
for j = 1:n
   % p = append(parameters.piper.md,list{j});
   p = [md, list{j}];
    A{21} = p;

   
    fid = fopen(parameters.piper.piperfileread, 'r+');
    
    for I = 1:24
        fprintf(fid, '%s\n', A{I});
    end
    fclose(fid);
   
   
        save('xyzDataMat','-v7.3');
        
    [Energies,nx,ny,nz]= testpiperread(matfri,parameters);
    EnergiesVec = [EnergiesVec Energies];
     clear Energies;
end
    matfri = k * parameters.piper.rotationwidth;
       EnergiesMat = [EnergiesMat;EnergiesVec]; 
       
    end
     

fname = append('EnergiesMat','_','dim',num2str(parameters.GRFModel.numofeigenvalues),'_xyz');
save (fname, 'EnergiesMat','-v7.3');
end
    