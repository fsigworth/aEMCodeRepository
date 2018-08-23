% function [fscs,vols,moi]=reGatherAndReconstruct(si,ri,iter,iTwin,gRoi,oldMoi,logs)
function roi=reGatherRois(ri,iTwin,gRoi,logs)
% Read the group roi files and combine them, writing out the composite roi.
% Do the reconstruction, write mrc volume files and the group moi files.
% Compute the fscs between our mrc volumes and the other twin's files.

groupFlags=reGetGroupFlags(ri,iTwin);
rois=cell(ri.nGroups,1);
rois{1}=gRoi;
if ri.nGroups>1
    mdisp(logs,['Gathering roi files ' datestr(now)]);
    for ig=2:ri.nGroups
        roiName=[reGetNameCode(ri,iter,iTwin,ig) 'roi.mat'];
        mdisp(logs,[' reading group: ' roiName datestr(now)]);
        [s,ok]=reCheckAndLoadMat(ri.tempPath, roiName,ri.timeout(2));
        if ~ok
            mdisp(logs,['Couldn''nt find the file ' ri.tempPath roiName datestr(now)]);
            error(['Couldn''nt find the file ' ri.tempPath roiName]);
        end;
        rois{ig}=s.gRoi;
    end;
    mdisp(logs,['all found.' datestr(now)]);
    
end;
roi=reCollectStructureFields(rois,groupFlags);

if ri.flags.saveRoi
    roiName=[ri.outPath reGetNameCode(ri,iter,iTwin) 'roi.mat'];
    mdisp(logs,['Saving the roi file ' roiName]);
    save(roiName,'roi');
end;

% reShowLatentVars(imgs,refs,ri,roi,iter,[7 10],['iter ' num2str(iter)]);
% drawnow;
%



function [moi,vols,fscs]=reGroupReconstruct(ri,iter,iTwin,iGroup,roi)
if iGroup>1
    roiName=[ri.outPath reGetNameCode(ri,iter,iTwin) 'roi.mat'];
    mdisp(logs,['Loading the roi file ' roiName]);
    load(roiName);
end;

nRefs=size(ri.angles,1);
rk=struct;  % make a fake ri structure for reGetGroupFlags
rk.nTwin=nRefs;
rk.nGroups=ri.nGroups;
gFlags=reGetGroupFlags(rk,1,iGroup);  % get the flags for splitting up the classes
ngRefs=sum(gFlags);
classNorms=roi.classNorms(:,:,gFlags,:);
classMeans=roi.classMeans(:,:,gFlags,:);
angles=ri.angles(gFlags,:);

mdisp(logs,'Reconstructions');
tic
n=ri.nCrop;
gProjs=zeros([n n ngRefs 2 ri.nVols],'single');

for iVol=1:ri.nVols
    realNorms=zeros(n,n,ngRefs,'single');
%     nCls=size(roi.classNorms,3);
    for i=1:ngRefs
        realNorms(:,:,i)=fftshift(real(ifftn(ifftshift(classNorms(:,:,i,iVol)))));
    end;
    %             classMeanSym=reSymmetrizeClasses(classMeans,ri);
    gProjs(:,:,:,1,iVol)=classMeans(:,:,:,iVol);
    gProjs(:,:,:,2,iVol)=realNorms;
end;
% usePar=ri.flags.useParFor %  && (ri.nVols>2);  % don't use parFor if few volumes
gFvo=reFourierInsertion(gProjs,angles,ri.symmetry,0);
            fvName=[reGetNameCode(ri,iter,iTwin,iGroup) 'fvo.mat'];
             mdisp(logs,[' saving group: ' fvName datestr(now)]);
             save(fvName,'gFvo');


if iGroup==1  % Collect all the Fourier volumes
    allFVols=gFvo;
    if ri.nGroups>1
        mdisp(logs,['Gathering fvols ' datestr(now)]);
        for ig=2:ri.nGroups
            fvName=[reGetNameCode(ri,iter,iTwin,ig) 'fvo.mat'];
            mdisp(logs,[' reading group: ' fvName datestr(now)]);
            [s,ok]=reCheckAndLoadMat(ri.tempPath, fvName,ri.timeout(2));
            if ~ok
                mdisp(logs,['Couldn''nt find the file ' ri.tempPath fvName datestr(now)]);
                error(['Couldn''nt find the file ' ri.tempPath fvName]);
            end;
            allFVols=allFVols+gFvo;
        end;
    mdisp(logs,['all found.' datestr(now)]);
    


%
%     Normalization of volumes
%     k=median(roi.varNs)/4;
k=mean(roi.imgAmps.^2)*ri.kFactor;  % mult by fsc/(1-fsc) says Sjors
mdisp(logs,'k',k);
otherVolNames=cell(ri.nVols,1);
vols=zeros(n,n,n,ri.nVols,'single');
moi.refVols=zeros(n,n,n,ri.nVols,'single');

for iVol=1:ri.nVols
    v=rsNormalizeReconstruction(allFVols(1,iVol),allFVols(2,iVol),k);
    
    vsd=sqrt(v(:)'*v(:)/numel(v));
    moi.refVols(:,:,:,iVol)=v*ri.volSD(iVol)*ri.refScale/vsd;
    vols(:,:,:,iVol)=v;
    
    volName=[reGetNameCode(ri,iter,iTwin,-iVol) '.mrc'];
    WriteMRC(vols(:,:,:,iVol),ri.pixA,[ri.outPath 'mrc/' volName]);
    otherVolNames{iVol}=[reGetNameCode(ri,iter,3-iTwin,-iVol) '.mrc'];    
end;
toc

%%
% Assign other model values
moi.sigmaN=sqrt(mean(roi.varNs));
moi.sigmaC=oldMoi.sigmaC;
moi.sigmaG=oldMoi.sigmaG;
moi.pRefs=oldMoi.pRefs;
moi.imgAmps=reSetImageAmplitudes2(si,ri,iTwin,roi.imgAmps);
moi.ringFits=oldMoi.ringFits+roi.ringFits;
moi.activeTrans=reAssignActiveTrans(roi.pTrans,ri.thresholds(1));
moi.pVols=mean(roi.pVols,2);
moi.intensiveFields={'sigmaN' 'sigmaC' 'sigmaG' 'refVols', 'pVols','pRefs'};

fracActiveTrans=sum(moi.activeTrans(:))/numel(moi.activeTrans);
mdisp(logs,'fracActiveTrans',fracActiveTrans);
%     fracActiveAlphas=sum(roi.pAlphas(:)>transThreshold)/numel(roi.pAlphas)
%     fracActiveRefs=sum(roi.pRefs(:)>transThreshold)/numel(roi.pRefs)

% Save the moi structure
if ri.nGroups>1 || ri.flags.saveMoi
    save([ri.outPath reGetNameCode(ri,iter+1,iTwin) 'moi.mat'],'moi');
end;

% Look up the other twin's volumes and compute the fsc
fscs=zeros(n/2-1,ri.nVols);
mdisp(logs,[datestr(now) ' Loading the other twin''s volumes ' otherVolNames{1}]);
for iVol=1:ri.nVols
    [altVol,s,ok]=reCheckAndLoadMRC([ri.outPath 'mrc/'], otherVolNames{iVol},ri.timeout(3),5);
    if ~ok
        mdisp(logs,[datestr(now) ' load failed']);
        break;
    else
        fscs(:,iVol)=FSCorr(altVol,vols(:,:,:,iVol));
    end;
end;

