function parameters = piperinitialization(parameters)

% Load piper dynamic shared library
% loadlibrary('./piper/libMatPiper')

% Initialize piper variables
piper.matfri = 0;
piper.matlri = 0;
%piper.output = int16(1);
piper.output = int16(0); % Flag false for do not output the data
%piper.piperfile = '../source/piper/piper_trypsin';

piper.Energies = 1;
piper.piperfile = 'piper_trypsin';

piper.DxFile = './piper/orig_grids/rec_grid.b4.l0.dx';
piper.md = '../../data/md_trypsin/';

parameters.piper = piper;
end
