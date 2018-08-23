function imgAmps=reSetImageAmplitudes2(si,ri,iTwin,roiImgAmps)
% Use the vesicle amplitudes to set image amplitudes: for each micrograph
% get the mean re-estimated image amplitude from roi.imgAmps; then
% scale the corresponding vesicle amplitudes to match.
twinInds=find(ri.twinFlags(:,iTwin));
mIndices=si.miIndex(twinInds);
ninds=numel(mIndices);
namps=numel(roiImgAmps);
if ninds ~= numel(roiImgAmps)
    error(['Inconsistent indices ' num2str([ninds numel(roiImgAmps)])]);
end;
imgAmps=zeros(namps,1,'single');
for i=1:numel(si.mi);  % loop over every micrograph
    p=mIndices==i;     % all the particles from it
    np=sum(p);
    if np>0
        miPAmps=roiImgAmps(p);
        avgPAmp=mean(miPAmps);
        miVAmps=si.sVesicle(p);
        normVAmp=mean(miVAmps);
        imgAmps(p)=miVAmps/normVAmp*avgPAmp;
    end;
end;
