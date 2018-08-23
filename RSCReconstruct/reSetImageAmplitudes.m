function imgAmps=reSetImageAmplitudes(si,roiImgAmps)
% Use the vesicle amplitudes to set image amplitudes: for each micrograph
% get the mean re-estimated image amplitude from roi.imgAmps; then
% scale the corresponding vesicle amplitudes to match.
nmi=numel(si.mi);
imgAmps=zeros(size(roiImgAmps),'single');
for i=1:nmi
    p=si.miIndex==i;
    np=sum(p);
    if np>0
        miPAmps=roiImgAmps(p);
        avgPAmp=mean(miPAmps);
        miVAmps=si.sVesicle(p);
        normVAmp=mean(miVAmps);
        imgAmps(p)=miVAmps/normVAmp*avgPAmp;
    end;
end;
