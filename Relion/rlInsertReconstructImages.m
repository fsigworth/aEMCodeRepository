% rlInsertReconstructImages.m
% Add the rlnReconstructImageName column to a subtracted particles star file.

% A typical line contains rlnImageName
%  000001@Refine3D/job217/Merged/001_X-1Y-1-2_v.mrcs 
% and the rlnReconstructImageName will be something like
%  000001@Extract/job220/Merged/001_X-1Y-1-2_u.mrcs 
% so we have to copy the image index, replace the path and the final character of the name.

newPath='Extract/job218/Merged/'
newSuffix='u'
nSuffChars=numel(newSuffix);

particleStarName='Refine3D/job214/run_data.star';
outStarName='Refine3D/job214/run_data_reconstr2.star';

[pnms,pdats]=ReadStarFile(particleStarName);

pdt=pdats{2};
np=numel(pdt.rlnImageName);
qdt=pdt;
qdt.rlnReconstructImageName=cell(np,1);
for i=1:np
    oldImgName=pdt.rlnImageName{i};
    atPtr=strfind(oldImgName,'@');
    [oldPath oldStack ext]=fileparts(oldImgName(atPtr+1:end));
    newStackName=[newPath oldStack(1:end-nSuffChars) newSuffix ext];
    qdt.rlnReconstructImageName{i}=[oldImgName(1:atPtr) newStackName];

end;

%%
qdats=[pdats(1); qdt];
disp(['Writing ' outStarName '...']);
WriteStarFile(pnms,qdats,outStarName);
disp('done.');