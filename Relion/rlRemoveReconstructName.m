% rlRemoveReconstructName.m
% Delete the reconstructImageName field from a particles file.
fName='rlnReconstructImageName';
prefix='su';

disp('Getting a particles.star file.');
[stName, stPath]=uigetfile('*.star','Select star file');
if isnumeric(stPath)  % user has clicked Cancel
    return
else
    cd(stPath);
end;
[stBlocks,stData,ok]=ReadStarFile(stName);

%%
st=stData{1};
if ~isfield(st,fName)
    disp(['No field found: ' fName]);
    return
end;
numParticles=size(st.rlnImageName,1);
disp([num2str(numParticles) ' particles.']);
st=rmfield(st,fName);
outName=[prefix stName];

disp(['writing ' outName '...']);
WriteStarFileStruct(st,'',outName);
disp('done.');
