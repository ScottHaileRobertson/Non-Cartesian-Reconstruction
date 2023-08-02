% To get this setup to run, you need a mex compiler:
%     In Linux:
%         edit mexopts.sh to have -sdt=c99
%         also line lstdc by http://stackoverflow.com/questions/17000903/mex-compiling-on-64-bit-linux-usr-bin-ld-cannot-find-lstdc

% Get current path
rootDistribDir = pwd();

% Add Recon package to path
disp('Adding Recon package to MATLAB path...');
path(genpath([rootDistribDir filesep() 'ReconPackage']),path);

% Compiling MEX code
disp('Compiling MEX code...');
Recon.SysModel.Proximity.mex_compile_proximity();
mex_compile_sparseGridder();
disp('Compilation worked');
