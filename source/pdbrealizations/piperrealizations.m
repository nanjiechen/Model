cd /projectnb/uqproj/StochasticDocking/Code/source/piper

loadlibrary('libMatPiper')

parameters.piper.matfri = 0;
parameters.piper.matlri = 0;
%piperfile = 'piper_argfile';
parameters.piper.piperfile = 'piper_trypsin';
parameters.piper.testnum = 2;
parameters.piper.Energies = 1;
parameters.piper.md = '/projectnb/uqproj/StochasticDocking/Code/data/md_trypsin/';
