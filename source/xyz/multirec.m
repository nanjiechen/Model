function [parameters] = multirec(parameters)
%  [status,cmdout] = system('ls ../data/md_trypsin/prot*.pdb');
% [status,cmdout] = system('ls ../data/md_trypsin/prot*.pdb > out');
 cd('../data/md_trypsin/');
 [status,cmdout] = system('ls prot*.pdb > out');
 [status,cmdout] = system('./minmax_list_pdb.pl out');
 [status,cmdout] = system('./minmax_list_pdb.pl out > multirecjoint.pdb');

% [status,cmdout] = system('ls ../data/PDBApprox/prot*.pdb');
% [status,cmdout] = system('ls ../data/PDBApprox/prot*.pdb > out2');
 cd('../PDBApprox/');
 [status,cmdout] = system('ls prot*.pdb > out2');
 [status,cmdout] = system('./minmax_list_pdb.pl out2');
 [status,cmdout] = system('./minmax_list_pdb.pl out2 > multirecjoint2.pdb');
 
[status,cmdout] = system('cat ../md_trypsin/multirecjoint.pdb multirecjoint2.pdb > multirec_trypsin_expand.pdb ');
[status,cmdout] = system('> multirecjoint2.pdb');
cd('../md_trypsin/');
[status,cmdout] = system('> multirecjoint.pdb');
 cd('../PDBApprox/');
[status,cmdout] = system('sed ''3d '' multirec_trypsin_expand.pdb >../../source/piper/multirec_trypsin_expand.pdb' );
 cd('/projectnb/uqproj/StochasticDocking/Code/source');
end