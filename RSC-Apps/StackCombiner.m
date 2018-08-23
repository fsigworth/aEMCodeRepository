% StackCombiner.m
% Merge multiple stacks into one.  We assume that pixA and image size are
% the same among the stacks.  A change in image size will give a fatal error.
% The output files are named like the last input file, but will have a '+',
% e.g.
% XXXX+tsi.mat
% XXXX+stall.mrc

% stackSuffix='stall.mrc';
stackSuffix='tstack.mrc';


% Get a list of tsi files to combine
[fname, pa]=uigetfile('*tsi.mat','Select tsi files','multiselect','on');
[rootPath, infoPath]=ParsePath(pa);
if ~iscell(fname)
    fname={fname};
end;
cd(pa);
%%
% Load the first one and use it as reference.
siName=fname{1};
disp(['Reading ' siName]);
load(siName);
nmis=numel(si.mi);
disp([num2str(nmis) ' micrographs']);
nims=numel(si.miIndex);  % number of particles
disp([num2str(nims) ' particles']);

p=strfind(siName,'tsi.mat');
stackName=[siName(1:p-1) 'stall.mrc'];
disp(['Reading ' stackName]);
stack=ReadMRC(stackName);
[n ny nstk]=size(stack);
disp(['Image size: ' num2str(n)]);
if nstk ~= nims
    error(['Inconsistent stack length: ' num2str([nstk nims])]);
end;

for index=2:numel(fname)
    sixName=fname{index};
    disp(['Reading ' sixName]);
    siTemp=load(sixName);
    siTemp=siTemp.si;
    nxmis=numel(siTemp.mi);
    disp([num2str(nxmis) ' micrographs']);

    %     Check for copies of ids, which means copies of micrographs
    ids=zeros(nmis,1);
    for i=1:nmis
        if numel(si.mi{i})>0  % not empty
            ids(i)=si.mi{i}.identifier;
        end;
    end;
    miOk=false(nxmis,1);
    for i=1:nxmis
        if numel(siTemp.mi{i})>0
            miOk(i)=~any(ids==siTemp.mi{i}.identifier);
        end;
    end;
    disp([num2str(sum(miOk)) ' micrographs to add']);
    nxims=numel(siTemp.miIndex);
    if nxims<1
        break;
    end;

%     Read the image stack
    p=strfind(sixName,'tsi.mat');
    stackxName=[sixName(1:p-1) 'stall.mrc'];
    disp(['Reading ' stackxName]);
    stackx=ReadMRC(stackxName);
    [n ny nxstk]=size(stackx);
    if nxstk ~= nxims
        error(['Inconsistent stack length: ' num2str([nstk nims])]);
    end;

    for j=1:nxmis  % loop through the micrographs
        if miOk(j)
            nmis=nmis+1;
            si.mi{nmis}=siTemp.mi{j};
            si.ctfs(:,:,nmis)=siTemp.ctfs(:,:,j);
            imIndices=find(siTemp.miIndex==j); % particles from this micrograph
            nxIndices=numel(imIndices);
            outIndices=(1:nxIndices)+nims;
            if nxIndices>0
                si.miIndex(outIndices)=nmis;
                si.miParticle(outIndices)=siTemp.miParticle(imIndices);
                si.alpha0(outIndices)=siTemp.alpha0(imIndices);
                si.yClick(outIndices)=siTemp.yClick(imIndices);
                si.rVesicle(outIndices)=siTemp.rVesicle(imIndices);
                si.sVesicle(outIndices)=siTemp.sVesicle(imIndices);
                stack(:,:,nims+1:nims+nxIndices)=stackx(:,:,imIndices);
                nims=nims+nxIndices
            end;
        end;
    end;
end;
si.ok=true(nims,1);

disp([num2str(nmis) ' Micrographs total']);
disp([num2str(nims) ' Images total']);

% Write out the composite stack
p=strfind(sixName,'tsi.mat');
if numel(p)>0
    outBasename=[sixName(1:p-1) '+'];
    siOutName=[outBasename 'tsi.mat'];
    stackOutName=[outBasename 'stall.mrc'];
    save(siOutName,'si');
    disp(['Writing ' siOutName]);
    WriteMRC(stack,si.pixA,stackOutName);
    disp(['Writing ' stackOutName]);
else
    disp('Nothing written');
end;