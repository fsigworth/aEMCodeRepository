function sVals=VesicleCircAmps(s,nTheta)
% function sVals=VesicleCircAmps(s,nTheta)

s=s(:);
% amplitude is a function of angle
    sExt=zeros(nTheta,1);
    ns=min(numel(s),nTheta);
    sExt(1:ns)=conj(s(1:ns));
    sVals=real(fft(sExt));
