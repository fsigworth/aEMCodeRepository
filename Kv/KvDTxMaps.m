% KvDTXMaps.m
% Try to make a difference map for Kv1.2 datasets with and without DTx

cd('/Users/fred/Documents/Documents - Katz/EMWorkOnDocs/Yangyu/KvDTX/Maps and model');
% [mTx,s]=ReadMRC('Kv1.2_DenTx_2.8A.mrc');
[mTx,s]=ReadMRC('Kv1.2DenTx_3.0A_bettertoxindensity.mrc');
mWt=ReadMRC('Kv1.2_WT_3.2A.mrc');
% mTx=SharpFilt(rsRotateImage(mTx,45),.3);
% mWt=SharpFilt(rsRotateImage(mWt,47),.3);

%% Extract the TM regions
    doWrite=0; % don't write out

n0=size(mTx,1);
ct0=ceil((n0+1)/2);
n1=144;
zsh=25;
zsh=27; % matches the model better.
tsh=1;
% Mask the T1 domain
mta=circshift(mTx,[0 0 zsh+tsh]);
t1MskZ=ones(n0,1);
t1MskZ(ct0+30:end)=0;
t1MskZ=GaussFilt(t1MskZ,.05);


t1Msk=shiftdim(repmat(t1MskZ,[1 n0 n0]),1);
% ShowSections(t1Msk.*mta)


mt0=Crop(t1Msk.*mta,n1);
mw=Crop(circshift(mWt,[0 0 zsh]),n1);

% Gotten from AlignSymmetricMaps:
% zsh1=-.102;
% rot1=4.2346;
zsh1=-.39; % values for 'bettertoxin' map.
rot1=2.74;

        fSh=FourierShift(n1,[0 0 zsh1]);
        mt1=real(ifftn(fftn(mt0).*fSh));
        mt=rsRotateImage(mt1,rot1);


figure(1);
ShowSections(SharpFilt(mt,.3),[],45);
set(gcf,'name','DenTx')
figure(2);
ShowSections(mw,[],45);
set(gcf,'name','WT');
%%
figure(3);
a=1.2;
ShowSections(a*SharpFilt(mt,.4)-mw,[],45)

%%
% fc0=.4;
% n3=numel(mt);
% F=ones(n3,2);
% 
% P=Simplex('init',fc0);
% for i=1:200
%     m1=SharpFilt(mt,P(1));
%     F(:,2)=m1(:);
%     as=LinLeastSquares(F,mw(:));
%     err=sum((F*as-mw(:)).^2);
%     imags(sum(as(2)*m1-mw,1));
%     title(num2str([i err P]));
%     drawnow;
% 
%     P=Simplex(err);
%     disp(err);
% end;
mwx=mw;
mtx=mt;
%% Find a linear filter (sum of Gaussians) to allow the best subtraction

ct=ceil((n1+1)/2);
t=6; % slab half-thickness
tm=12; % slab for projection image of masked maps
% Make a mask for the fitting region
mski=fuzzymask(n1,3,23,2,[ct ct 50]);
msko=fuzzymask(n1,3,35,4,[ct ct ct+10]);
msk=max(0,msko-mski);
% 2nd try
msko=fuzzymask(n1,3,[21 21 1000],2); % cylinder around pore.
mski=fuzzymask(n1,3,22,2,[ct ct 40]); % toxin blank
mskp=fuzzymask(n1,3,[4 4 1000],1); % cylinder around ions
mskp2=fuzzymask(n1,3,[7 7 50,2,[ct ct 90]]);
msk=msko.*(1-mski).*(1-mskp).*(1-mskp2);
% ShowSections(msk.*mtx);


mt=mtx;
mtm=msk.*mt;

% mw=SharpFilt(mwx,.33);
mw=mwx;
mwm=msk.*mw;
% fc values for the many parallel Gaussian filters.
fcs=[.6 .4 .37 .36 .35 .34 .33 .32 .3 .29 .28 .27 .26 .25 .24 .23 .22 .21 .2];
nTerms=numel(fcs);
n3=numel(mtm);
F=zeros(n3,nTerms+1);
F(:,end)=.05;
F1=F;
for i=1:nTerms
    fc=fcs(i);
    F(:,i)=reshape(SharpFilt(mtm,fc),n3,1);
    F1(:,i)=reshape(SharpFilt(mt,fc),n3,1);
end;
%-------Least squares fitting here------
a=lsqr(F,mwm(:));
disp(a);
mtmf=reshape(F*a,[n1 n1 n1]);
mtf=reshape(F1*a,[n1 n1 n1]);

figure(3);
mscl=1000;
mysubplot(231);
% imaga(sum(mtmf,1)*mscl+100);
imaga(mean(mtmf(ct-tm:ct+tm,:,:),1)*mscl+100);
title('Tx');
mysubplot(232);
imaga(mean(mwm(ct-tm:ct+tm,:,:),1)*mscl+100);
title('Wt');
mysubplot(233);
%     imaga(sum(mtmf-mwm,1)*mscl+100);
    imaga(mean(mtmf(ct-t:ct+t,:,:)-mwm(ct-t:ct+t,:,:),1)*mscl+100);
    title('Tx-Wt')
mysubplot(234)
% imaga(sum(mtf,1)*mscll+100);
imaga(mean(mtf(ct-t:ct+t,:,:),1)*mscl+100);
mysubplot(235);
imaga(mean(mw(ct-t:ct+t,:,:),1)*mscl+100);

mysubplot(236);
%     imaga(sum(mtf-mw,1)*mscl+100);
    imaga(mean(mtf(ct-t:ct+t,:,:)-mw(ct-t:ct+t,:,:),1)*mscl+100);

%% Write out the diff map with various amounts of filtering
pixA=s.pixA;
diffMap=mtf-mw;

if doWrite
diffMap33=SharpFilt(diffMap,.3*pixA);
diffMap4=SharpFilt(diffMap,.25*pixA);
diffMap5=SharpFilt(diffMap,.2*pixA);
WTMap=mw;
outPath='DiffMaps/';
        disp(['Writing files to ' outPath])
WriteMRC(diffMap,pixA,[outPath 'diffMap.mrc']);
WriteMRC(diffMap33,pixA,[outPath 'diffMap33.mrc']);
WriteMRC(diffMap4,pixA,[outPath 'diffMap4.mrc']);
WriteMRC(diffMap5,pixA,[outPath 'diffMap5.mrc']);
WriteMRC(WTMap,pixA,[outPath 'WTMap.mrc']);
WriteMRC(mt,pixA,[outPath 'matchTxMap.mrc'] );
WriteMRC(mtf,pixA,[outPath 'matchFiltTxMap.mrc'] );
disp('done.')
    end;




return


%% plot of slices through the maps.

% msko=fuzzymask(n1,3,[21 21 1000],2); % cylinder around pore.
% mski=fuzzymask(n1,3,22,2,[ct ct 40]); % toxin blank
% mskp=fuzzymask(n1,3,[4 4 1000],1); % cylinder around ions
% mskp2=fuzzymask(n1,3,[7 7 50,2,[ct ct 90]]);
% msk=msko.*(1-mski).*(1-mskp).*(1-mskp2);

figure(4);

scz=65;
scy=73;
for i=1:10; scz=95-5*i; subplot(10,1,i); plot([mtf(:,scy,scz) mw(:,scy,scz)]); ylabel(scz); end;

figure(5);
mtf2=mtf.*msk;
mw2=mw.*msk;
for i=1:10; scz=95-5*i; subplot(10,1,i); plot([mtf2(:,scy,scz) mw2(:,scy,scz)]); ylabel(scz); end;

figure(6);
ShowSections(mtf2);

return

%% make a mask based on the model.
[comp,prot]=SolventAndProteinDensity('WT_fit_WTmap_molrep_ABCD_R1.pdb');
%
protScl=DownsampleGeneral(circshift(prot,[0 0 0]),144,1/s.pixA);
% align the maps
% [rot,dz,protScla]=AlignSymmetricMaps(protScl,500*SharpFilt(mw,.2));
% turns out that rot=0.3 degrees, dz=-.1: negligible.
protScla=protScl;
protMsk=GaussFilt(GaussFilt(protScla,.1)>.1,.2);
ShowSections(protMsk.*diffMap);
% WriteMRC(protMsk.*diffMap,s.pixA,'MaskedDiffMap.mrc');
%% Make special masks to accomodate the DTx
txSphere=fuzzymask(144,3,17,1,[0 0 -27]+73);
txMsk=min(1,protMsk+txSphere);
ShowSections(SharpFilt(diffMap,.25).*txMsk);
WriteMRC(txSphere.*diffMap,s.pixA,'MaskedTx.mrc');
WriteMRC(txMsk.*diffMap,s.pixA,'MaskedDiffMap.mrc');
WriteMRC(mw.*protMsk,s.pixA,'MaskedWT.mrc');
WriteMRC(mt.*txMsk,s.pixA,'MaskedKvTx.mrc');
