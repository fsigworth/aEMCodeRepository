function [si,imgs]=reMakeFakeData(sim,mapName)
    n=sim.n;
    vesRs=[300 300];
    if mod(sim.pIso,1)==0  % if all are iso or rso, simulate only that.
        isos=sim.pIso;
    else
        isos=[0 1];
    end;
    [ri, angles0]=reSetRefAngles(sim.angleSteps,sim.angleLimits,isos,false);
    na=size(angles0,1);
    nim=na*sim.nMi;
    angles=repmat(angles0,1,1,sim.nMi);
    angles=reshape(permute(angles,[3 1 2]),nim,3);

    si=struct;
    si.isos=rand(nim,1)<sim.pIso;
    angles(:,1)=180*si.isos+sim.sigmaA*randn(nim,1);
    si.pixA=4.988;
    si.mbnOffset=9.62;
    si.sVesicle=.008*ones(nim,1,'single');
    si.miIndex=ones(nim,1,'uint16');
    si.miParticle=(1:nim)';
    si.alpha0=single(360*rand(nim,1)-180);
    si.rVesicle=single(vesRs(1)+(vesRs(2)-vesRs(1))*rand(nim,1));
    effR=si.rVesicle+sign(si.isos-.5)*si.mbnOffset;
    sins=sind(angles(:,2));
    si.yClick=sins.*effR;
    si.shifts=sim.sigmaC*randn(nim,2);
    si.weights=[1 1];
    si.ctfs=ones(n,n,1,'single');
    si.ctfs=abs(CTF(n,si.pixA,.025,2,2,0,.02));
    si.sVesicle=0.01*ones(nim,1,'single');
    si.mi=struct;
    si.activeFlags=true(nim,1);
    si.activeFlagLog={[date '  reMakeFakeData']};
    %
    % ctfsc=zeros(nc,nc,nMicrographs,'single');
    % for j=1:nMicrographs
    %     d=sm.defoci(j);
    %     si.ctfs(:,:,j)=abs(CTF(n,pixA,lambda,d,2,B0+B1*d,ctAlpha));
    %     ctfsc(:,:,j)=abs(CTF(nc,pixA,lambda,d,2,B0+B1*d,ctAlpha));
    % end;
    % sm.isos=rand(nImgs,1)<pIso;
    origVol=arGetRefVolumes(si.pixA,n,mapName,1);
    noise=sim.sigmaN*randn(n,n,nim,1,1);
    signal0=sim.imgScale*reMakeTemplates(origVol,angles);  % refs(x,y,iRef,iVol,iTwin)
    fShifts=FourierShift([n n],si.shifts);
    signalf=real(ifft2(fft2(signal0).*fShifts.*repmat(ifftshift(si.ctfs),1,1,nim)));
    imgs=signalf+noise;

