% rlExtendGroupClusters.m
% Expand groups to local clusters of micrographs, chosen to yield a minimum
% set of particles.

minClusterParticles=200;
% Load a particles.star file from classification or selection.
% 
inStarPath='RSC9/Class3D295_sel/';
inStarName='particles_iso_all.star';
% pStarName='Class3D/job187/run_it025_data.star';

outStarPath='RSC9/Class3D295_regroup/';
CheckAndMakeDir(outStarPath,1);
outStarName=inStarName;

% We'll find the entries iVes in the ves-part file corresponding to each line of
% the particles.star file.
%%
inName=[AddSlash(inStarPath) inStarName];
disp(['Reading ' inName]);
[pnm,pdat]=ReadStarFile(inName);
opt=pdat{1};
pts0=pdat{2};
%%
nParticles=numel(pts0.rlnMicrographName);

% Get the sorted group names
[gpNames,~,gpParticles]=unique(pts0.rlnGroupName);
% operate on the basis of group names, assign group numbers.
nGroups=numel(gpNames);
gpIndices=1:nGroups;

% Get sorted micrograph names
[micNames,nameInds,partInds]=unique(pts0.rlnMicrographName);
nMics=numel(micNames);
disp([num2str(nMics) ' unique micrographs']);
%
pts=pts0;

groupIndexIncrement=0;
lastImic=0;
nMicGroups=zeros(nMics,nGroups);
szGroups=0;
disp('Original group  addedGroups');
for igrp=1:nGroups
    grpParticles=pts0.rlnGroupNumber==igrp;
    newGroupParticles=0;
    for iMic=1:nMics % look at each particle micrograph name
        % Want to find partInds
        mgParticles=find(grpParticles & partInds==iMic);
%        mgParticles=find(pts.rlnGroupNumber(partInds==iMic)==igrp);
        nmgp=numel(mgParticles);
        nMicGroups(iMic,igrp)=nmgp;
        if nmgp>0
            newGroupParticles=newGroupParticles+nmgp;
%                          disp([iMic nmgp newGroupParticles]);

            pts.rlnGroupNumber(mgParticles)=igrp+groupIndexIncrement;
            gpName=pts.rlnGroupName{nameInds(iMic)};
            newGroupName=[gpName '+' num2str(groupIndexIncrement)];
            pts.rlnGroupName(mgParticles)={newGroupName};
            if newGroupParticles>=minClusterParticles
                lastIncrement=groupIndexIncrement;
                lastGroupName=newGroupName;
                lastGroupParticles=newGroupParticles;
                szGroups(igrp+groupIndexIncrement,1)=newGroupParticles;

                groupIndexIncrement=groupIndexIncrement+1;
                newGroupParticles=0;
                lastImic=iMic;
            end;
        end;
    end;
    if newGroupParticles<minClusterParticles % combine with previous group
%         display(['Backing up at group ' num2str(igrp) ' to micrograph ' num2str(lastImic)]);
%         display(['newGroupParticles ' num2str(newGroupParticles) ' nmgp ' num2str(nmgp)]);
        groupIndexIncrement=lastIncrement;
        newGroupParticles=lastGroupParticles;
        for iMic=lastImic+1:nMics
            mgParticles=grpParticles & partInds==iMic;
            nmgp=sum(mgParticles);
            newGroupParticles=newGroupParticles+nmgp;
%             disp([iMic nmgp newGroupParticles])
            if nmgp>0
                pts.rlnGroupNumber(mgParticles)=igrp+lastIncrement;
                pts.rlnGroupName(mgParticles)={lastGroupName};
            end;
        end;
%                                 disp([iMic nmgp newGroupParticles]);

    end;
    disp([igrp groupIndexIncrement]);
end;
%
        hdrText=['# version 30001, ' num2str(igrp+groupIndexIncrement) ' groups. from ' inName]
        fullOutName=[outStarPath outStarName];
        disp(['Writing ' fullOutName]);
         WriteStarFile(pnm,{opt pts},fullOutName,hdrText);

