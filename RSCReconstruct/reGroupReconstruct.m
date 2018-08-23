function [moi,vols]=reGroupReconstruct(ri,si,iter,iTwin,roi,moi,gFvs,logs)
% function [moi,vols]=reGroupReconstruct(ri,si,iter,iTwin,roi,oldMoi,gFvs,logs)
% Call this on the iGroup==1 job.  Gathers and sums the Fourier volumes, and makes
% the reconstruction set vols.  Also assigns the model structure moi.
tCyc=1;  % cycle time

allFVols=gFvs;
if ri.nGroups>1
    mdisp(logs,[datestr(now) ' Gathering fvols']);
    for ig=2:ri.nGroups
        fvName=[reGetNameCode(ri,iter,iTwin,ig) 'fvs.mat'];
        mdisp(logs,[datestr(now) ' reading group: ' fvName]);
        [gFvs,ok]=reCheckAndLoadMat(ri.tempPath, fvName,ri.timeout(2),tCyc,1);
        if ~ok
            mdisp(logs,['Couldn''nt find the file ' ri.tempPath fvName datestr(now)]);
            error(['Couldn''nt find the file ' ri.tempPath fvName]);
        end;
        for i=1:numel(gFvs)
            allFVols(i).PadFT=gFvs(i).PadFT+allFVols(i).PadFT;
        end;
    end;
    mdisp(logs,[datestr(now) ' all found.']);
end;

%
%     Normalization of volumes
%     k=median(roi.varNs)/4;
mdisp(logs,'Volume normalization');

k=mean(roi.imgAmps.^2)*ri.kFactor;  % mult by fsc/(1-fsc) says Sjors
if iter<ri.highResIter
    k=k*4;
end;

mdisp(logs,'k',k);

n=gFvs.n;
vols=zeros(n,n,n,ri.nVols,'single');  % returned volume
for iVol=1:ri.nVols
    v=rsNormalizeReconstruction(allFVols(1,iVol),allFVols(2,iVol),k);
    vols(:,:,:,iVol)=v;
end;
moi.refVols=reNormalizeModels(ri,vols);
%%
% % Assign other model values
% moi.sigmaN=sqrt(mean(roi.varNs));
% moi.sigmaC=oldMoi.sigmaC;
% moi.sigmaG=oldMoi.sigmaG;
% moi.pRefs=oldMoi.pRefs;
% moi.imgAmps=reSetImageAmplitudes2(si,ri,iTwin,roi.imgAmps);
% moi.ringFits=oldMoi.ringFits+roi.ringFits;
% %moi.activeTrans=reAssignActiveTrans(roi.pTrans,ri.thresholds(1));
% [moi.activeTrans,moi.activeAlphas,moi.activeRVs]=reAssignActivePars(ri,roi);
% moi.pVols=mean(roi.pVols,2);
% moi.intensiveFields={'sigmaN' 'sigmaC' 'sigmaG' 'refVols', 'pVols','pRefs'};
% 
% fracActiveTrans=sum(moi.activeTrans(:))/numel(moi.activeTrans);
% fracActiveAlphas=sum(moi.activeAlphas(:))/numel(moi.activeAlphas);
% fracActiveRVs=sum(moi.activeRVs(:))/numel(moi.activeRVs);
% mdisp(logs,'fracActiveTrans',fracActiveTrans);
% mdisp(logs,'fracActiveAlphas',fracActiveAlphas);
% mdisp(logs,'fracActiveRVfs',fracActiveRVs);
% %     fracActiveAlphas=sum(roi.pAlphas(:)>transThreshold)/numel(roi.pAlphas)
% %     fracActiveRefs=sum(roi.pRefs(:)>transThreshold)/numel(roi.pRefs)


