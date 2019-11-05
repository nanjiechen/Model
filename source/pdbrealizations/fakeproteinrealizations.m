% Create fake realizations for purposes of debugging

function [parameters,data] = fakeproteinrealizations(methods,parameters);

numofrealizations = 3;
numofgridpoints = 10;
mu = 0;
sig = 1;
data.protein = normrnd(mu,sig,[numofgridpoints,numofrealizations]);


