function moi=reEstimateModelValues(ri,si,roi,iTwin,oldMoi,allActive,logs)

moi=oldMoi;

% Assign other model values
moi.sigmaN=sqrt(mean(roi.varNs));
moi.sigmaC=oldMoi.sigmaC;
moi.sigmaG=oldMoi.sigmaG;
moi.pRefs=oldMoi.pRefs;
moi.imgAmps=reSetImageAmplitudes2(si,ri,iTwin,roi.imgAmps);

if allActive
    ri.thresholds=[0 0 0];
end;
[moi.activeTrans,moi.activeAlphas,moi.activeRVs]=reAssignActivePars(ri,roi);
moi.pVols=mean(roi.pVols,2);
moi.intensiveFields={'sigmaN' 'sigmaC' 'sigmaG' 'refVols', 'pVols','pRefs'};

fracActiveTrans=sum(moi.activeTrans(:))/numel(moi.activeTrans);
fracActiveAlphas=sum(moi.activeAlphas(:))/numel(moi.activeAlphas);
fracActiveRVs=sum(moi.activeRVs(:))/numel(moi.activeRVs);
mdisp(logs,'fracActiveTrans',fracActiveTrans);
mdisp(logs,'fracActiveAlphas',fracActiveAlphas);
mdisp(logs,'fracActiveRVfs',fracActiveRVs);
%     fracActiveAlphas=sum(roi.pAlphas(:)>transThreshold)/numel(roi.pAlphas)
%     fracActiveRefs=sum(roi.pRefs(:)>transThreshold)/numel(roi.pRefs)
