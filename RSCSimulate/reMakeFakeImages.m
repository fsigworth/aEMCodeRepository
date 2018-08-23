function [si,imgs]=reMakeFakeImages(sim,mis, mapName)
% Create noiseless projection images.  The sim structure has fields
% n,nmi,npm      image size, number of micrographs, particles per micrograph
% mis, mbnOffsetA  cell array of mi structures; membrane offset in A
% angleLimits  e.g. [10 20 90] for |alpha|<10, beta 20...160, gamma 0...90
% pIso, sigmaC    iso probability; click sigma in ds pixels
% rVes, sVes      min, max values for radius (in ds pixels), amplitude
% aParticle       min, max of particle scaling relative to vesicle
% mergeMode       merging mode: 1=normal; 2=ctf of first only; 3=ctf of
% first, and unflipped.
% We use a simple effCTF calculation, but include a typical PW filter.

usePWFilter=sim.usePWFilter;

    n=sim.n;
    nim=sim.nmi*sim.npm;

    si=struct;
    [angles,isos]=rsGetRandomAngles(nim,sim.angleLimits,sim.pIso);
    si.pixA=sim.ds*mis{1}.pixA;
    si.mbnOffset=sim.mbnOffsetA/si.pixA;
    
    si.sVesicle=(sim.sVes(1)+(sim.sVes(2)-sim.sVes(1))*rand(nim,1));
    si.miIndex=reshape(repmat(uint16(1:sim.nmi),sim.npm,1),nim,1);
    si.miParticle=repmat(uint16(1:sim.npm)',sim.nmi,1);
    si.alpha0=single(360*rand(nim,1)-180);
    si.rVesicle=single(sim.rVes(1)+(sim.rVes(2)-sim.rVes(1))*rand(nim,1));
    effR=si.rVesicle+sign(isos-.5)*si.mbnOffset;
    sines=sind(angles(:,2));
    si.yClick=sines.*effR;
    si.shifts=sim.sigmaC*randn(nim,2);
    si.ctfs=ones(n,n,sim.nmi,'single');
    for i=1:sim.nmi
        if usePWFilter
            pwf=meGetNoiseWhiteningFilter(mis{i},n,sim.ds,1,.001*si.pixA);
        else
            pwf=1;
        end;
        si.ctfs(:,:,i)=pwf.*meGetEffectiveCTF(mis{i},n,sim.ds,sim.mergeMode);
    end;
    si.sVesicle=single(sim.sVes(1)+(sim.sVes(2)-sim.sVes(1))*rand(nim,1));
    si.mi=mis;
    si.weights=mis{1}.weights;
    si.activeFlags=true(nim,1);
    si.activeFlagLog={[date '  reMakeFakeData']};
    %
    si.origVol=arGetRefVolumes(si.pixA,n,mapName,1);
    disp(['Making ' num2str(nim) ' projections...']);
    signal0=reMakeTemplates(si.origVol,angles);  % refs(x,y,iRef,iVol,iTwin)
    disp('Filtering...');
    fShifts=FourierShift([n n],si.shifts);
    ctfs=zeros(n,n,nim,'single');
    si.sParticle=zeros(nim,1,'single');
    si.angles=angles;
    si.isos=isos;
    aParticle=sim.aParticle(1)+(sim.aParticle(2)-sim.aParticle(1))*rand(nim,1);
    for i=1:nim
        s=si.sVesicle(i)*aParticle(i);
        si.sParticle(i)=s;
        ctfs(:,:,i)=s*ifftshift(si.ctfs(:,:,si.miIndex(i)));
    end;
    imgs=real(ifft2(fft2(signal0).*fShifts.*ctfs));
    disp('done.');

