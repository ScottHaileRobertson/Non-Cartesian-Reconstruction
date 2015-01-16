% Get current path
rootDistribDir = pwd();

% Load 3p libs
disp('Loading 3rd party libs...');
path(genpath([rootDistribDir filesep() '3pLibs' filesep() 'AutoLoad']),path);

% Add Recon package to path
disp('Adding Recon package to MATLAB path...')
path(genpath([rootDistribDir filesep() 'ReconPackage']),path);

% Compiling MEX code
disp('Compiling MEX code...')

