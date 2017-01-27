% Get the starting directory
rootDistribDir = pwd();

% Change to the mex directory
cd([rootDistribDir filesep() 'ReconPackage' ]);

% Compile mex code
mex -largeArrayDims mex_thread_calcSparseGridMatrix.c
mex -largeArrayDims mex_thread_sparseMultiply.c;

% Change back to starting directory 
cd(rootDistribDir);