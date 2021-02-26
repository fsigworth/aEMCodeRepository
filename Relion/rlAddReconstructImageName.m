% rlAddReconstructImageName
% add the column to a particles.star file

inStarName='Select/job191/particles.star';
outStarName='Select/job191/particles+unsub.star';

disp(['Reading ' inStarName]);
[nm,dat]=ReadStarFile(inStarName);
imgNames=dat{2}.rlnImageName;

whos imgNames
%%
outNames=imgNames;
nim=numel(imgNames);
for i=1:nim
    [partInd,micName]=rlDecodeImageName(imgNames{i});
    [pa,nm,ex]=fileparts(micName);
    nm(end)='u'; % swap the final character from 'v' to 'u'
    name=imgNames{i};
    if strndcmp(name(end-4:end),'.mrcs')
        name(end-5)='u'; % mark unsubtracted
        disp(name);
        return
    end;
end;