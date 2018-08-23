function [si,imgs]=rsStackMerge(sia,imga)
% Merge the array of stackInfo structures sia into the
% accumulated structure si. This undoes the rsStackSplit operation.
% Images can also be merged if the cell array of
% image stacks imga is also given.  We do no checking for inconsistent
% entries, but depend on the origParticle and origMi fields for
% identification.
% 
%
ns=numel(sia);
nim=0;
nmi=0;
for i=1:ns
    nim=max(nim,max(sia(i).origParticle));
    nmi=max(nmi,max(sia(i).origMi));
end;

si=struct;
si.origParticle=zeros(nim,1,'uint32');  % zero where nothing is filled in.

si.origMi=cell(1,nmi);

for i=1:ns
    
%     per-particle information
    pinds=sia(i).origParticle;
si.origParticle(pinds)=pinds;
si.miIndex(pinds)=sia(i).origMi(sia(i).miIndex);
si.miParticle(pinds,i)=sia(i).miParticle;
si.alpha0(pinds,1)=sia(i).alpha0;
si.yClick(pinds,1)=sia(i).yClick;
si.rVesicle(pinds,1)=sia(i).rVesicle;
si.sVesicle(pinds,1)=sia(i).sVesicle;
if isfield(sia,'ctfGroupIndex')
    si.ctfGroupIndex(pinds,1)=sia(i).ctfGroupIndex;
end;

% per-mi information
    minds=sia(i).origMi;
    si.origMi=zeros(nim,1,'uint16');
    si.origMi(minds)=minds;
    si.mi(minds)=sia(i).mi;
    si.ctfs(:,:,minds)=sia(i).ctfs;    

end;

activeFlags=si.origParticle>0;

% global
si.pixA=sia(1).pixA;
si.weights=sia(1).weights;
if isfield(sia,'mbnOffset')
    si.mbnOffset=sia(1).mbnOffset;
end;
if isfield(sia,'ctfGroups')
    si.ctfGroups=sia(1).ctfGroups;
end;

% images
imgs=[];
if nargin>1 && nargout>1 % we are going to operate on images
    nx=size(imga{1},1);
    imgs=zeros(nx,nx,nim,'single');
    for i=1:ns
        imgs(:,:,sia(i).origParticle)=imga{i};
    end;
        % Now that we've created a sparse structure, we shrink it thusly:
    [si,imgs]=rsStackSplit(activeFlags,si,imgs);
else
    si=rsStackSplit(activeFlags,si);
end;

