%% Main file to run, allows for multiple scenarios to be run in series
% Requires download and install the following: 
% 
% Distmesh 
%    https://popersson.github.io/distmesh/index.html
% CVX 
%    http://cvxr.com/cvx/download/
% MEaSUREs matlab plug in + data files
%   https://www.mathworks.com/matlabcentral/fileexchange/47329-measures
% BedMachine matlab plug in + data files
%   https://www.mathworks.com/matlabcentral/fileexchange/69159-bedmachine


clear
close all
clc

addpath('lib') 
addpath('grid') 

if(~ismac && ~ispc)
    % Start-up business on sherlock is hard. You must have a MATLAB folder and 
    % startup.m script to get everything set up. Duplicate this from your local 
    addpath('lib') 
    addpath ../MATLAB/
    run startup.m
end

% Name of scenarios to run, only 1 map file used here.
nameToRun = ["ISSM"];
mapsToRun = ["strainMesh035.mat"];
config    = defaultConfig();

for j = 1:length(mapsToRun)
    for i = 1:length(nameToRun)
        clearvars -except nameToRun mapsToRun i j config
        str = nameToRun(i);
        mapFile = mapsToRun(j);
        modelRunner;
    end
end
