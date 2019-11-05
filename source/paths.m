d = fileparts(pwd);
a = genpath(d);
path(path,a);

% Remove folders from paths where piper is built and tested

rmpath('../data/triset/')

rmpath('../toolbox/libmol2/')

rmpath('../toolbox/libgrid2/')

rmpath('../toolbox/librestraints/')

rmpath('../toolbox/piper1/')
