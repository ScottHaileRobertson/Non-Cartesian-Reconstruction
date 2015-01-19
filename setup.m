% Get current path
rootDistribDir = pwd();

% Add Recon package to path
disp('Adding Recon package to MATLAB path...');
path(genpath([rootDistribDir filesep() 'ReconPackage']),path);

% Compiling MEX code
disp('Compiling MEX code...');
Recon.SysModel.Proximity.mex_compile_proximity();
