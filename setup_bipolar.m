% code to setup bipolar environment and data paths

currentPath = pwd;
parentPath = fileparts(currentPath);
%setenv('BIPOLAR_DATA', fullfile(currentPath, 'data'));
setenv('BIPOLAR_DATA', '/scratch/bipolar_expedition');
addpath(genpath(fullfile(currentPath)))
