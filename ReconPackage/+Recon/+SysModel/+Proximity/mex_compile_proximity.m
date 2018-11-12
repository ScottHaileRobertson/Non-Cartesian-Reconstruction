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
        mex -n  -largeArrayDims sparse_gridding_distance_mex.c;
        db_inplace(mfilename,'COMPILER PROBLEMS, We ran in test mode, use terminal commands preceding this line as examples to compile manually. Of note, -ansi may cause trouble on linux systems.');
    end
end
% Change back to starting directory
cd(rootDistribDir);