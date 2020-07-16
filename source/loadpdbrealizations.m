%% Load all the realizations

pdbfilenames = dir('../data/cleanset/*.pdb')

numofrealizations = size(pdbfilenames,1);

for n = 1 : numofrealizations,
    
    tic;
    temp = pdbread(pdbfilenames(n).name);
 
    x = vertcat(temp.Model.Atom.X );  
    X(1:length(x),n) = x;
    
    y = vertcat(temp.Model.Atom.Y );  
    Y(1:length(y),n) = y;
    
    z = vertcat(temp.Model.Atom.Z );  
    Z(1:length(z),n) = z;
    toc;
        
end