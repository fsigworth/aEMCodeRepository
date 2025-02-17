% HRPickingModel.m
% Create 3D models with and without the micelle
% of Kv TM region and the complete alpha4
% for high-res picking.

targetPixA=1.09;
computeCCs=0;
doWrite=0; % Write out the models?
shiftToMatchModel=0;
n1=144; % working map box size
% n2=256; % Output map size, matches input map
ct1=ceil((n1+1)/2);
% cd('/Users/fred/EMWork/Yangyu/20210224_YW_sel'); % on Katz
cd('/Users/fred/Documents/Documents - Katz/EMWorkOnDocs/Yangyu/20210224_YW_sel/'); % on Mini2
% cd('/Volumes/EMWork/Yangyu/20210224_YW_sel'); % Katz mounted
pdbName='W366F_v4.2(delete_loop add ion)_VSDchecked.pdb'; % TM only
outPath='/Volumes/D255/20181216/No5Graphene/sq05_1/HRRefs/';

% [mp0,s]=ReadMRC('postprocess_masked_2.55A.mrc');
% mCtr=CenterOfMass(mp0)
% mp=circshift(mp0,round(-mCtr));

[mp0,s]=ReadMRC('Postprocess/job171/postprocess.mrc');
n2=size(mp0,1);
ct2=ceil((n2+1)/2);
baseRotDegrees=0;
% mpr=rsRotateImage(mp,baseRotDegrees);
mp1=DownsampleGeneral(mp0,n1,s.pixA); % convert to 1A/pixel
figure(1);
% ShowSections(mp1,[102 82 73]);
ShowSections(mp1,[],45);
drawnow;

%% Make density from the model
[dComp,dProt]=SolventAndProteinDensity(pdbName);
dComp=Crop(dComp-dComp(1),n1);
dProt=Crop(dProt,n1);
n1=size(dComp,1);

dSolv=dComp-dProt;

%% Shift the exp map to match the model, which has only TM region.
figure(3);
cc3=fftshift(real(ifftn(fftn(mp1).*conj(fftn(dProt)))));
ShowSections(cc3);
[val,zShift]=max3di(cc3) % shift is 21 pixels in z
mp1sh=circshift(mp1,round(ct1-zShift)); % Center the tm region
figure(1);
ShowSections(mp1sh);
%% Calculate the relative amplitude of the exp map.
amp=(dProt(:)'* mp1sh(:))/(dProt(:)'*dProt(:))
% value is .0035

%% Get micelle 1D density
ySlab=squeeze(mean(mp1sh(:,ct1-12:ct1+9,:),2));
xySlab=mean(ySlab(ct1+32:ct1+38,:),1);
figure(2);
mysubplot(121);
imags(ySlab);
mysubplot(122);
plot(xySlab);

% micDens=circshift(xySlab(:),[21 0]);
micDens=xySlab(:);
md0=micDens(ct1+34);
micDens(1:ct1-31)=md0; % 73-31
micDens(ct1+31:end)=md0; % 73+31
% Symmetrize the profile
micDens(ct1:-1:2)=micDens(ct1:end);
micDens1=GaussFilt(micDens,.15);
plot(micDens1);
%% Replace protein with modeled micelle density to get the micelle map

protMask=GaussFilt(GaussFilt(dProt,.05)>.01,.2);
micelleDens=repmat(shiftdim(micDens1,-2),[n1 n1 1]); % spread the micelle density
figure(3);
ShowSections(GaussFilt( (1-protMask).*mp1sh+protMask.*micelleDens,.05),[ct1 ct1 ct1+28]);

% Create a mask that clips everything above and below the TM region.
t1Mask=zeros(n1,n1,n1,'single');
w=26;
t1Mask(:,:,ct1-w:ct1+w)= 1;

t1Maskf=GaussFilt(t1Mask,.05);
ShowSections(GaussFilt( (1-protMask).*mp1sh+protMask.*micelleDens,.2),[ct1 ct1 ct1+28]);
drawnow;
% Create the micelle model density. Outside the protein, use mp1sh.
% Inside the protein, use micelleDens.
mMap0=GaussFilt( (1-protMask).*mp1sh.*t1Mask+protMask.*micelleDens*1,.1);
ShowSections(mMap0,[ct1 ct1 ct1+28]);
% drawnow;
% ShowSections(t1Maskf.*mp1sh);

drawnow;
%

% Get the subtracted TM region alone
tMapsh=(t1Maskf.*mp1sh-mMap0).*protMask;
ShowSections(tMapsh,[],45);
% Get the subtracted alpha-subunits alone
aMapsh=mp1sh-mMap0;



%% Unshift the maps
t1Msk1=circshift(protMask.*t1Maskf,round(-zShift-ct1));
tMap1=circshift(tMapsh,round(zShift-ct1)); % undo the shift of the TM region
ShowSections(tMap1,[],45);
aMap1=circshift(aMapsh,round(zShift-ct1)); % whole alpha subunit, micelle subtracted
micMap1=circshift(mMap0,round(zShift-ct1)); % undo the micelle shift too.

% Try to solvent-flatten the alpha-subunit map
protMsk1=circshift(protMask,round(zShift-ct1));
aMsk1=GaussFilt(max(GaussFilt(aMap1,.05)>.006,protMsk1),.05); % Mask that includes T1 domains
amMap1=aMap1-1.2*GaussFilt(aMap1,.004); % undo most of the overshoot
% aMod1=GaussFilt(aMap1(:,:,[133 136]),.03); % get extra density at the top

% Identify the T1 region
isT1=GaussFilt(GaussFilt((GaussFilt(aMap1,.05)>.006) > protMsk1,.1)>.8,.1);

compMap1=amMap1.*min((protMsk1+isT1),1); % overshoot-removed, masked alpha subunit.
compMap1shx=circshift(Crop(compMap1,n2),round(ct1-zShift));

compWeakMap1=t1Msk1.*max(0,2*(aMap1-.012))+tMap1;
compWeakMap1shx=circshift(Crop(compWeakMap1,n2),round(ct1-zShift));


if doWrite % save the 1A pixel maps, all centered on the TM region: micelle, TM, whole alpha
    save HRPicking/Figs/TM_Micelle_maps.mat mMap0 tMapsh compMap1shx micDens1
end;


tMap2=DownsampleGeneral(tMap1,n2,1/targetPixA); % TM. Convert back to the orig pixel size
% tmaMap2=DownsampleGeneral(amMap1,n2,1/targetPixA); % TM plus T1 domain: whole alpha subunit
compMap2=DownsampleGeneral(compMap1,n2,1/targetPixA); % TM+T1
compWeakMap2=DownsampleGeneral(compWeakMap1,n2,1/targetPixA); % TM+T1
mMap2=DownsampleGeneral(micMap1,n2,1/targetPixA); % micelle

zShift2=CenterOfMass(tMap2); % get the CM of just the TM region
tMap2u=circshift(tMap2,round(-zShift2));
compMap2u=circshift(compMap2,round(-zShift2));

%%
tMag=tMap2(:)'*tMap2(:)

% aMap2=DownsampleGeneral(aMap1,n2,1/s.pixA);
% amMap2=DownsampleGeneral(amMap1.*aMsk1,n2,1/s.pixA);
figure(4);
% ShowSections(amMap2,[ct2 ct2 ct2+10],45);
ShowSections(compMap2,[ct2 ct2 ct2+10],45);

compMag=compMap2(:)'*compMap2(:)
% Correct the micelle map also.
mMap2=DownsampleGeneral(micMap1,n2,1/s.pixA);
mtMag=tMap2(:)'*mMap2(:)

% Show correlation of TM with micelle, and ACF of TM, and ACF of alpha
mtCC=fftshift(real(ifftn(fftn(mMap2).*conj(fftn(tMap2)))));
tmAC=fftshift(real(ifftn(abs(fftn(tMap2)).^2)));
compAC=fftshift(real(ifftn(abs(fftn(compMap2)).^2)));

figure(5);
ShowSections(mtCC); % This is the cross-correlation of micelle with protein. The
% peak is about 3 while the autocorrelation of the protein peak is 75. The
% alpha peak is about 130.

figure(6);
ShowSections(tmAC);
mysubplot(3,3,1);
title('TM region autocorrelation');

figure(7);
ShowSections(compAC);
mysubplot(3,3,1);
title('alpha composite ACF');

%% ---------Write the output files---------
if doWrite
    CheckAndMakeDir(outPath);
    disp(['In ' outPath ' :']);
    WriteMRC(tMap2,targetPixA,[outPath 'tmMap.mrc']); % Map of TM Region
    disp('  tmMap.mrc');
    WriteMRC(mMap2,targetPixA,[outPath 'micMap.mrc']); % Micelle map
    disp('  micMap.mrc');
    WriteMRC(compMap2,targetPixA,[outPath 'compMap.mrc']);
    disp('  compMap.mrc');
    WriteMRC(compWeakMap2,targetPixA,[outPath 'compWeakMap.mrc']);
    disp('  compWeakMap.mrc');
    disp('...written.');
else
    disp('Nothing writtten.');
end;

return
%%
if computeCCs
    tMag=tMap2(:)'*tMap2(:)
    
    % aMap2=DownsampleGeneral(aMap1,n2,1/s.pixA);
    % amMap2=DownsampleGeneral(amMap1.*aMsk1,n2,1/s.pixA);
    figure(4);
    % ShowSections(amMap2,[ct2 ct2 ct2+10],45);
    ShowSections(compMap2,[ct2 ct2 ct2+10],45);
    
    compMag=compMap2(:)'*compMap2(:)
    % Correct the micelle map also.
    mtMag=tMap2(:)'*mMap2(:)
    
    % Show correlation of TM with micelle, and ACF of TM, and ACF of alpha
    mtCC=fftshift(real(ifftn(fftn(mMap2).*conj(fftn(tMap2)))));
    tmAC=fftshift(real(ifftn(abs(fftn(tMap2)).^2)));
    compAC=fftshift(real(ifftn(abs(fftn(compMap2)).^2)));
    
    figure(5);
    ShowSections(mtCC); % This is the cross-correlation of micelle with protein. The
    % peak is about 3 while the autocorrelation of the protein peak is 75. The
    % alpha peak is about 130.
    
    figure(6);
    ShowSections(tmAC);
    mysubplot(3,3,1);
    title('TM region autocorrelation');
    
    figure(7);
    ShowSections(compAC);
    mysubplot(3,3,1);
    title('alpha composite ACF');
    
end;

return
%% Which is the best postprocess map?
% mp=circshift(mp,round(mCtr));


d1=dir('Postprocess');
for i=3:numel(d1)
    if d1(i).isdir
        nm=['Postprocess/' d1(i).name '/postprocess_masked.mrc']
        mp2=ReadMRC(nm);
        diff=mp2(:)-mp0(:);
        diff2=sqrt(diff'*diff)
    end;
end;
