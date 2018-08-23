function [mAccum nAccum, mAccumOdd, nAccumOdd]=deGetAccumSum(folderPath, skipFrames)
% Given a directory containing DE movie folders, compute the raw means of each
% movie, store it in the directory, and return the global accumulation.
% If the directory only contains
% ACCUM Also return the accumulation
% of every other movie.
fp=AddSlash(folderPath);
mAccum=0;
nAccum=0;
mAccumOdd=0;
nAccumOdd=0;
accumName=[fp 'AccumSum.mat'];
accumOddName=[fp 'AccumSumOdd.mat'];
odd=1;
r=FindFilenames(fp,'DE',1);  % Get dir names starting with DE
for i=1:numel(r)
    movieDir=[fp AddSlash(r{i})];
    [m nim]=deGetRawSum(movieDir,skipFrames);
    if nim>0
        save([movieDir 'isum.mat'],'m','nim');
    end;
    mAccum=mAccum+m;
    nAccum=nAccum+nim;
    if odd
        mAccumOdd=mAccumOdd+m;
        nAccumOdd=nAccumOdd+nim;
    end;
    odd=1-odd;
end;
if nAccum>0
    save(accumName,'mAccum','nAccum');
    save(accumOddName,'mAccumOdd','nAccumOdd');
end;

