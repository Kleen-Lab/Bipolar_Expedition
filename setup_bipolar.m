% code to setup bipolar environment and data paths

currentPath = pwd;
parentPath = fileparts(currentPath);
setenv('BIPOLAR_DATA', '/scratch/bipolar_expedition'); % modify for user 
setenv('RESULTS','/scratch/bipolar_expedition/results'); % modify for u ser 
addpath(genpath(fullfile(currentPath)))
