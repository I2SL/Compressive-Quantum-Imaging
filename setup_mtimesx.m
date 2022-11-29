% use if downloading functions from scratch
try
    mtimesx('SPEEDOMP',['OMP_SET_NUM_THREADS(',num2str(maxNumCompThreads('automatic')),')']);
catch
    if exist('mtimesx','dir') ~=7
        unzip('https://github.com/nonkreon/mtimesx/archive/master.zip','mtimesx');
        addpath(genpath([pwd filesep 'mtimesx' filesep 'mtimesx-master' filesep 'src']));
        savepath;
    end
    try
        mtimesx('SPEEDOMP',['OMP_SET_NUM_THREADS(',num2str(maxNumCompThreads('automatic')),')']);
    catch
        mex -setup;
        mtimesx('SPEEDOMP',['OMP_SET_NUM_THREADS(',num2str(maxNumCompThreads('automatic')),')']);
    end
end