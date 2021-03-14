% rlRestoreParticleMetadata

% Read a particles.star file A (such as one we made)
%  and another one B having selected particles.
% Make C, a copy of B, except
%  - copy the optics group header from A to C
%  - copy metadata from A for each particle in B (recognized by same micrograph
% and coordinates) and place this into C. The fields to be replaced are
theFields={'rlnOpticsGroup'
           'rlnGroupName'
           'rlnGroupNumber'};

starAName='RSC/particleAll9_intens+frac_7505.star';
starBName='Select/job236/particles.star';
starNewBName='Select/job236/particles_new_metadata.star';

%%
disp(['Reading ' starAName]);
[nmA,datA]=ReadStarFile(starAName);
ptsA=datA{2};
nAParticles=numel(ptsA.rlnMicrographName);
%%
    disp(['Reading ' starBName]);
        [nmB,datB]=ReadStarFile(starBName);
        bPts=vdat{2};
    end;
%% Copy B into C
nmC=nmB;
datC=datB;
    %% Copy the optics info from A
    datC{1}=datA{1};

%%  ----needs work
[bNames,~,bPartInds]=unique(bPts.rlnMicrographName);
bnp=numel(bNames);
disp([num2str(bnp) ' unique micrographs']);

[aNmes,~,aPartInds]=unique(aPts.rlnMicrographName);
disp([num2str(numel(aNames)) ' micrograph refs']);
pvMap=zeros(bnp,1);
iVes=zeros(nParticles,1);
for i=1:bnp % look at each particle micrograph name
    % find all the particles having that name
    iParts=find(partInds==i);
    pxs=pts.rlnCoordinateX(iParts);
    pys=pts.rlnCoordinateY(iParts);
    % find the correspoinding vesicle name index
    j=find(strcmp(partNames{i},vesNames),1);
    vesRange=find(vesInds==j);
%     if numel(vesRange)>1
%         break;
%     end;
    % now we find which coordinates match
    vxs=ves.ptlX(vesRange);
    vys=ves.ptlY(vesRange);
    dists=hypot(pxs-vxs',pys-vys');
    [minDists,iVesLocal]=min(dists,[],2);
    iVes(iParts)=vesRange(iVesLocal);
 end;

 % At this point, iVes(i) gives the line in ves cooresponding to pts(i)
 %%
 totalRso=0;
 nClasses=max(pts.rlnClassNumber);
 for icls=1:nClasses
 class=icls;
 
 sel=pts.rlnClassNumber==class;
 
 psis=pts.rlnAnglePsi(sel);
 psis=mod(psis,360);
 vesPsi=-ves.vesPsi(iVes(sel));
 vesPsi=mod(vesPsi,360);
 
 vesRs=ves.vesR(iVes(sel));

 figure(3);
 plot(vesPsi,psis,'.');
 xlabel('Vesicle psi');
 ylabel('Particle psi');
 title(['Class ' num2str(class)]);
 
 figure(4);
 psiDiff=mod(vesPsi-psis+90,360);
 hist(psiDiff,100);
 iso=psiDiff<180;
 rso=psiDiff>180; 
 rsoFraction(icls)=sum(rso)/(sum(rso)+sum(iso));
 title(['Class ' num2str(icls) ' rso-fraction ' num2str(rsoFraction(icls))]);
 disp([icls rsoFraction(icls)]);
 totalRso=totalRso+sum(rso);
pause(1);
 end;
 disp(pStarName);
totalRso
 return
 
 
 %% Create an edited star file
 selectRso=0;
 outStarName='RSC/iso_from_Class3D_121_data.star';
 psiResidual=pts.rlnAnglePsi+ves.vesPsi(iVes);
 psiDiff=mod(90-psiResidual,360);
 rso=psiDiff>180;
 hist(psiDiff,100);
 title(['Total ' num2str(sum(rso)) ' RSO Particles of ' num2str(numel(rso))]);
 disp(['Writing ' num2str(sum(rso)) ' particles to ' outStarName '...']);
 hdrText='# version 30001 from rlMatchVesicleAngles.m';
 if selectRso
  pFlags={[] rso}; 
    disp('Selecting RSO...');
 else
     pFlags={[] ~rso};
    disp('Selecting ~RSO...');
 end;
 
  WriteStarFile(pnm,pdat,outStarName,hdrText,pFlags);
disp('done.'); 
 
 