% Get the starting directory
rootDistribDir = pwd();

% Change to the mex directory
cd(fullfile(rootDistribDir , 'ReconPackage' ));

do_compile=1;
if exist('mex_thread_calcSparseGridMatrix.mexa64','file')...
    && exist('mex_thread_sparseMultiply.mexa64','file')
    f1=dir('mex_thread_calcSparseGridMatrix.c');
    f2=dir('mex_thread_sparseMultiply.c');
    f1o=dir('mex_thread_calcSparseGridMatrix.mexa64');
    f2o=dir('mex_thread_sparseMultiply.mexa64');
    if f1.datenum < f1o.datenum ...
            && f2.datenum < f2o.datenum
        do_compile=0;
    end
end
if do_compile
    try
        % Compile mex code
        mex -largeArrayDims mex_thread_calcSparseGridMatrix.c
        mex -largeArrayDims mex_thread_sparseMultiply.c;
    catch err
        warning(err.message);
        mex -n -largeArrayDims mex_thread_calcSparseGridMatrix.c
        mex -n -largeArrayDims mex_thread_sparseMultiply.c;
        db_inplace(mfilename,'COMPILER PROBLEMS, We ran in test mode, use terminal commands preceding this line as examples to compile manually. Of note, -ansi may cause trouble on linux systems.');
    end
end

% Change back to starting directory 
cd(rootDistribDir);