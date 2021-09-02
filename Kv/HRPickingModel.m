% HRPickingModel.m
% Create 3D model and 2D projections for high-res picking.

n1=144; % working map box size
n2=196; % Output map size
ct1=ceil((n1+1)/2);
cd('/Users/fred/EMWork/Yangyu/20210224_YW_sel');
pdbName='W366F_v4.2(delete_loop add ion)_VSDchecked.pdb'; % TM only

% [mp0,s]=ReadMRC('postprocess_masked_2.55A.mrc');
% mCtr=CenterOfMass(mp0)
% mp=circshift(mp0,round(-mCtr));

[mp0,s]=ReadMRC('Postprocess/job171/postprocess.mrc');
baseRotDegrees=0;
% mpr=rsRotateImage(mp,baseRotDegrees);
mp1=DownsampleGeneral(mp0,n1,s.pixA); % convert to 1A/pixel
figure(1);
% ShowSections(mp1,[102 82 73]);
ShowSections(mp1,[],45);
drawnow;

%%
[dComp,dProt]=SolventAndProteinDensity(pdbName);
dComp=Crop(dComp-dComp(1),n1);
dProt=Crop(dProt,n1);
n1=size(dComp,1);

dSolv=dComp-dProt;

%%
figure(3);
cc3=fftshift(real(ifftn(fftn(mp1).*conj(fftn(dProt)))));
ShowSections(cc3);
[val,zShift]=max3di(cc3) % shift is 21 pixels in z
mp1=circshift(mp1,round(ct1-zShift)); % Center the tm region
figure(1);
ShowSections(mp1);
%%
amp=(dProt(:)'* mp1(:))/(dProt(:)'*dProt(:))
% value is .0035

%% Get micelle 1D density
ySlab=squeeze(mean(mp1(:,ct1-12:ct1+9,:),2));
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
%%
protMask=GaussFilt(GaussFilt(dProt,.05)>.01,.2);
micelleDens=repmat(shiftdim(micDens1,-2),[n1 n1 1]);
figure(3);
ShowSections(GaussFilt( (1-protMask).*mp1+protMask.*micelleDens,.05),[ct1 ct1 ct1+28]);

t1Mask=zeros(n1,n1,n1,'single');
w=26;
t1Mask(:,:,ct1-w:ct1+w)= 1;

t1Maskf=GaussFilt(t1Mask,.05);
ShowSections(GaussFilt( (1-protMask).*mp1+protMask.*micelleDens,.2),[ct1 ct1 ct1+28]);
drawnow;
mMap0=GaussFilt( (1-protMask).*mp1.*t1Mask+protMask.*micelleDens*1,.1);
ShowSections(mMap0,[ct1 ct1 ct1+28]);
% drawnow;
% ShowSections(t1Maskf.*mp1);

drawnow;
%
ShowSections(t1Maskf.*mp1-mMap0,[],45);
%
tMap0=(t1Maskf.*mp1-mMap0).*protMask;
ShowSections(tMap0,[],45);
%%
tMap1=circshift(tMap0,round(zShift-ct1));
ShowSections(tMap1,[],45);

mMap1=circshift(mMap0,round(zShift-ct1));

tMap2=DownsampleGeneral(tMap1,n2,1/s.pixA);
tmMag=tMap2(:)'*tMap2(:)
mMap2=DownsampleGeneral(mMap1,n2,1/s.pixA);
mtMag=tMap2(:)'*mMap2(:)
mtCC=fftshift(real(ifftn(fftn(mMap2).*conj(fftn(tMap2)))));
tmAC=fftshift(real(ifftn(abs(fftn(tMap2)).^2)));
figure(5);
ShowSections(mtCC); % This is the cross-correlation of micelle with protein. The
% peak is about 3 while the autocorrelation of the protein peak is 75.
figure(6);
ShowSections(tmAC);
mysubplot(3,3,1);
title('TM region autocorrelation');

%% Write the output files
outPath='HRPicking/';
CheckAndMakeDir(outPath);
WriteMRC(tMap2,s.pixA,[outPath 'tmMap.mrc']);
WriteMRC(mMap2,s.pixA,[outPath 'micMap.mrc']);


%%
% See the total power in acf as a function of resolution.
freqs2=(0:n2/2-1)'/(n2*s.pixA);
sp2=RadialPowerSpectrum(tMap2);
cumPower=cumsum(freqs2.*sp2);
sps=sp2;
cumPowers=cumPower;
bVals=[0 50 100]; % first one is assumed to be zero.
labels={'B=0'};
for k=2:numel(bVals)
B=bVals(k);
gaussB=exp(-freqs2.^2*B/2);
sps(:,k)=sp2.*gaussB;
cumPowers(:,k)=cumsum(freqs2.*sp2.*gaussB);
labels{k,1}=['B=' num2str(B)];
end;
mysubplot(3,3,8);
semilogy(freqs2,sps);
legend(labels);
title('Spectral density');

mysubplot(3,3,9);
plot(freqs2,cumPowers);
legend(labels);
title('Image power');
% text(0.35,cumPower(find(freqs2>.35,1)),['B=' num2str(B)],'verticalalignment','bottom');
% hold off;
% our final picking reference is tMap2. The final micelle model is tMap2


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