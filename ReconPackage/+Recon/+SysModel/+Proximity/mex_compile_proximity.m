% Gvet the starting directory
rootDistribDir = pwd();

% Change to the mex directory
cd(fullfile(rootDistribDir , 'ReconPackage' ,...
    '+Recon' , '+SysModel' , '+Proximity' ));

% Compile mex code
do_compile=1;
if exist('sparse_gridding_distance_mex.mexa64','file')
    f1=dir('sparse_gridding_distance_mex.c');
    f2=dir('sparse_gridding_distance.c');
    fo=dir('sparse_gridding_distance_mex.mexa64');
    if f1.datenum < fo.datenum ...
            && f2.datenum < fo.datenum
        do_compile=0;
    end
end
if do_compile
    try
        mex -largeArrayDims sparse_gridding_distance_mex.c;
    catch err
        warning(err.message);
        fprintf('cd %s\n',pwd);
        mex -n  -largeArrayDims sparse_gridding_distance_mex.c;
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