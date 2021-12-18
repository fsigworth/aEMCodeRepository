% rlMatchVesicleAngles.m

loadVesData=0;
insertVesiclePsis=0;

% Load a particles.star file, may contsain a selected subset from Refine3D etc.
% pStarName='Refine3D/job140/run_data.star';
pStarName='Refine3D/job056/run_data.star';  % 2o0211122 dataset

% Load a vesicle-particle file from the entire dataset.
vStarName='RSC1/ves_particles_v.star';


% We'll find the entries iVes in the ves-part file corresponding to each line of
% the particles.star file.
%%
disp(['Reading ' pStarName]);
[pnm,pdat]=ReadStarFile(pStarName);
pts=pdat{2};
nParticles=numel(pts.rlnMicrographName);
%%
if loadVesData || ~exist('ves','var')
    disp(['Reading ' vStarName]);
    if strndcmp(vStarName,'.mat',4)
        load(vStarName); % loads ves structure
    else
        [vnm,vdat]=ReadStarFile(vStarName);
        ves=vdat{2};
    end;
end;

nVesicles=numel(ves.vesMicrographName);
%%
[partNames,~,partInds]=unique(pts.rlnMicrographName);
npn=numel(partNames);
disp([num2str(npn) ' unique micrographs']);

[vesNames,~,vesInds]=unique(ves.vesMicrographName);
disp([num2str(numel(vesNames)) ' vesicle refs']);
pvMap=zeros(npn,1);
iVes=zeros(nParticles,1);
disp('Matching particle locations...');
for i=1:npn % look at each particle micrograph name
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
disp('done.');

 % At this point, iVes(i) gives the line in ves cooresponding to pts(i)
 %%
 bins=0:4:359; % histo bins
 nBins=numel(bins);
     nClasses=max(pts.rlnClassNumber);
 h=zeros(nBins,nClasses);
 fStr=sprintf('%%%uu',floor(log10(nParticles)+1)); % format string
 for icls=1:nClasses
     class=icls;
     clsString=['Class ' num2str(icls)];
     if nClasses>1
         sel=pts.rlnClassNumber==class;
     else
         sel=true(size(pts.rlnMicrographName)); % select everything
     end;

     psis=pts.rlnAnglePsi(sel);
%      psis=mod(psis,360);
     vesPsi=ves.vesPsi(iVes(sel));
%      vesPsi=mod(vesPsi,360);
     psiDiff=mod(-vesPsi-psis+90,360);
     h(:,icls)=hist(psiDiff,bins);
     iso=psiDiff<180;
     rso=psiDiff>180;
     fprintf(['%s ' fStr ' iso ' fStr ' rso, fraction %.3f \n'],clsString,sum(iso),sum(rso),sum(rso)/numel(rso));
 end;

      psis=pts.rlnAnglePsi;
     vesPsi=ves.vesPsi(iVes);
     
    if insertVesiclePsis
        disp('Inserting the vesiclePsi field.');
     pts.vesiclePsi=vesPsi;
     pdat{2}=pts;
    end;

    psiDiff=mod(-vesPsi-psis+90,360);
     iso=psiDiff<180;
     rso=psiDiff>180;
     fprintf(['\n Total: ' fStr ' iso ' fStr ' rso, fraction %.3f \n\n'],sum(iso),sum(rso),sum(rso)/numel(rso));

     figure(4);
     bar(bins,h,'stacked');
     legend;
     ylabel('Frequency')
     xlabel('Psi angle difference');
     title([num2str(sum(iso)) ' ISO     ' num2str(sum(rso)) ' RSO']);
 
     [outPath]=AddSlash(fileparts(vStarName));
     ok=MyInput('Write a star file? ',0);
     while ok
         selRSO=MyInput('RSO (-1 for both orientations) ? ',1);
         if selRSO<0
             orString='riso'
             orFlags=true(nParticles,1);
         elseif selRSO
             orString='rso'
             orFlags=rso;
         else
             orString='iso'
             orFlags=iso;
         end;
         classSel=MyInput('Classes or 0?',0)
         if all(classSel)==0
             clsString='all';
             flags=true(nParticles,1) & orFlags;
         else
             clsString=sprintf('%u',classSel)
             flags=any(pts.rlnClassNumber==classSel,2) & orFlags;
         end;
         numParticles=sum(flags);
         path=input(['output file path [''' outPath ''']?']);
         if numel(path)>0
             outPath=AddSlash(path);
             CheckAndMakeDir(outPath,1);
         end;
         outStarName=['particles_' orString '_' clsString '.star'];
         outName=input(['File name? ''' outStarName ''' ? '])
         if numel(outName)<5
             outName=outStarName
         end;
        hdrText=['# version 30001, ' orString ' particles, class ' clsString ' from ' pStarName]
        totalParticles=sum(flags)
        fullOutName=[outPath outName];
        disp(['Writing ' fullOutName]);
         WriteStarFile(pnm,pdat,fullOutName,hdrText,{[] flags});
         
         ok=MyInput('Write another star file? ',0);
     end;
