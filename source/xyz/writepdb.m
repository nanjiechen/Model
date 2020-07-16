
%Obtain Energies
function [parameters] = writepdb(parameters)
l = parameters.stochasticmodel.numofgrids;
d = '/projectnb/uqproj/StochasticDocking/Code/data/md_trypsin/';
d2 = '/projectnb/uqproj/StochasticDocking/Code/data/PDBApprox/'; 
% PDBStruct = pdbread('prot_001_out_nmin.pdb');
fid = fopen('prot_001_out_nmin.pdb','r');
C = textscan(fid,'%s %s %s %s %s %s %s %s %s %s %s');  
input.recordName = {C{1}{4:end-1}};
input.atomName = {C{3}{4:end-1}};
input.resName = {C{4}{4:end-1}};
input.chainID = {C{5}{4:end-1}};
resNum = {C{6}{4:end-1}};
resNum = sprintf('%s ', resNum{:});
input.resNum = sscanf(resNum, '%f');
occupancy = {C{10}{4:end-1}};
occupancy = sprintf('%s ', occupancy{:});
input.occupancy = sscanf(occupancy, '%f');
betaFactor = {C{11}{4:end-1}};
betaFactor = sprintf('%s ', betaFactor{:});
input.betaFactor= sscanf(betaFactor, '%f');

x = parameters.stochasticmodel.xRecMap;
y = parameters.stochasticmodel.yRecMap;
z = parameters.stochasticmodel.zRecMap;
n = size(x,1);
for i = 1 : l
% sx = num2cell(x(:,i)');
% sy = num2cell(y(:,i)');
% sz = num2cell(z(:,i)');
% [PDBStruct.Model.Atom.X]= sx{:};
% [PDBStruct.Model.Atom.Y]= sy{:};
% [PDBStruct.Model.Atom.Z]= sz{:};
input.X = x(:,i);
input.Y = y(:,i);
input.Z = z(:,i);
File = append(d2,'prot_',num2str(i.','%04d'),'_approximation_nmin.pdb');
parameters.piper.approxPDB = File;
input.outfile = File;
mat2pdb(input);
end
fclose(fid);
end