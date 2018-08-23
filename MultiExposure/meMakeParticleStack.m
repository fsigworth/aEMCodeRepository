% function [stackInfo stack ctfs]=meMakeParticleStack(miArray,ds,nb,indexLimits)
displayOn=1;
subtractVesicles=1;
ds=2;
nb=48;  % box size
if displayOn
    figure(1);
    SetGrayscale;
end;

% definition of stackInfo structure
%     si.miBaseFilename(miIndex)
%     si.pixA         % Copied from runInfo
%     si.miIndex(i)   % Index of the micrograph corresponding to particle i (index of the micrographInfos)
%     si.miParticle(i)% Index of the particle within the micrograph.  This allows lookup of x,y,vesR etc.
%     si.alpha0(i)    % alpha rotation applied to the image
%     si.yClick(i)    % Copied from the mi struct.  xClick is zero after rotation.
%     si.rVesicle(i)  % Copied from the mi struct
basePath='/Volumes/TetraData/EMWork/Hideki/121118/PC1PC2liposome_SUP_14k/';
infoPath='Info/';
stackPostfix='stack48.mat';

miArray=meLoadInfoFiles([basePath infoPath]);

nmi=numel(miArray);
disp([basePath infoPath]);
disp([num2str(nmi) ' files found.']);

npar=0;
pixA0=miArray{1}.pixA;
for j=1:nmi  % Scan first to count the particles
    mi=miArray{j};
    if pixA0 ~= mi.pixA
        error(['Inconsistent pixA value at ' mi.baseFilename 'mi.mat']);
    end;
    ok=mi.particle.type>0;
    npar=npar+sum(ok);
end;
disp([num2str(npar) ' particles total.']);
pixA=ds*pixA0
si.pixA=pixA;
si.miBaseFilename=cell(nmi,1);
si.miIndex=single(zeros(npar,1));
si.miParticle=single(zeros(npar,1));
si.alpha0=single(zeros(npar,1));
si.yClick=single(zeros(npar,1));
si.rVesicle=single(zeros(npar,1));

stack=single(zeros(nb,nb,npar));
ctfs=single(zeros(nb,nb,nmi));

nst=0;  % number of particles scanned
fs=RadiusNorm(nb)/pixA;
for j=1:nmi
    mi=miArray{j};
    %     Read the merged image, subtract vesicles and apply PW filter.
    imageFile=[basePath mi.procPath mi.baseFilename 'm.mrc'];
    disp(imageFile);
    m=ReadEMFile(imageFile);
    if subtractVesicles
        ms=m-meMakeModelVesicles(mi,size(m),0,1); % use CTF
    else
        ms=m;
    end;
    H=meGetNoiseWhiteningFilter(mi,size(m));
    mf=real(ifftn(fftn(ms).*ifftshift(H)));
    if displayOn
        subplot(1,2,1);
        imacs(GaussFilt(m,.1));
        title(mi.baseFilename,'interpreter','none');
        subplot(1,2,2);
        imacs(GaussFilt(mf,.1));
        drawnow;
    end;
    %     Get the CTF and PW filter corresponding to a particle image
    ct1=meGetEffectiveCTF(mi,nb,ds);
    ct2=meGetNoiseWhiteningFilter(mi,nb,ds);
    ok=mi.particle.type>0;
    indices=find(ok);
    nok=numel(indices);
    %     Assign the stack, ctf and si fields
    stack(:,:,nst+1:nst+nok)=meGetParticleStack(mf,mi,nb,indices);
    ctfs(:,:,j)=ct1.*ct2;
    si.miBaseFilename{j}=mi.baseFilename;
    si.miIndex(nst+1:nst+nok)=j;
    si.miParticle(nst+1:nst+nok)=indices(:);
    disp([num2str(nok) ' particles']);
    nst=nst+nok;
end;
%% Write out the stack data after the last merged image
save([basePath mi.procPath mi.baseFilename stackPostfix], 'si','stack','ctfs');