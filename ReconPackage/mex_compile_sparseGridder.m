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
        wrds= sprintf('%s\n\t%s\n\t%s\n\t%s\n\t%s\n\t%s',... 
            'COMPILER PROBLEMS, We ran in test mode!', ...1
            'Use terminal commands preceding this line', ...2
            'as examples to compile manually. ', ...3
            'Of note, -ansi may cause trouble on linux systemms.', ...4
            'DONT FORGET TO chage to correct directory ', ...5
            pwd ... 6
            );
        db_inplace(mfilename,wrds);
    end
end

% Change back to starting directory 
cd(rootDistribDir);