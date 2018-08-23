% TestPickingForRSC3

% 
TestFakeData;

% Make eigenCCs
sigmaN=1;
s0=.06;  % image scale factor
imgCTF=.06*real(ifftn(fftn(img).*ifftshift(abs(CTF(n,1,.025,3,2,300,.07)))));
m=imgCTF+sigmaN*randn(n,n);
m=m+fuzzymask(n,2,10,3,[75 129]);
figure(1);
subplot(1,1,1);
imacs(m);
ccrefs=zeros(n,n,nterms);

for i=1:nterms
    ccrefs(:,:,i)=ifftshift(Crop(eigenims(:,:,i),n));
end;

ccs=zeros(n,n,nterms);
fimg=fftn(m);
disp('Computing ccs');
for i=1:nterms
    ccs(:,:,i)=real(ifftn(fimg.*conj(fftn(ccrefs(:,:,i)))));
end;

ptrs=rsGetAngleMap(n,r0-membraneOffset,inds);
%%
disp('Assembling image');
% want to make an nfactors x n x n x nGamma array of vectors
pArray=ptrs(1:n,1:n);
vArray=shiftdim(vlist(:,pArray,:),2);  % nGamma * nterms * n^2
ccArray=shiftdim(reshape(ccs,n*n,nterms),1);
projamps=reshape(projamps,nAngles,nGamma*nMirror+nrds);
paArray=reshape(projamps(pArray,1:nGamma),n*n,nGamma);
mPower=sum(ccs.^2,3);  % noise part will be of magnitude sigmaN^2

cimg=zeros(nGamma,n*n);
for i=1:n*n
    cimg(:,i)=vArray(:,:,i)*ccArray(:,i);
end;

% norm=repmat(mPower(:)',nGamma,1);
nccs=cimg'./repmat(sqrt(mPower(:)),1,nGamma);
ncc=reshape(max(nccs,[],2),n,n);

cimgn=cimg./(s0*paArray');  % expected amplitudes = 1.  scc
[mxv mxi]=max(cimgn);  % normalized by expected value.

mximg=reshape(mxv,n,n);
mxind=reshape(mxi,n,n);


figure(2);
imacs(mximg>0.7);
figure(3);
pimg=reshape(paArray,n,n,nGamma);
plot([sect(mximg) sect(ncc) sect(mxind-1)*.05]);
%  plot([sect(mximg) sect(ncc)]);
disp('done.')

blankRadius=nt*.5;
maxAmp=1.5;
minAmp=0.6;

mximg2=mximg;
nfound=0;
[amp ix iy]=max2di(mximg2);
amps=[];
coords=[];
figure(2);
while amp>minAmp
    % blank the cross-correlation
    mximg2=mximg2.*(1-fuzzymask(n,2,blankRadius,0.1,[ix iy]));
    if amp <= maxAmp && amp > minAmp  % valid peak
        nfound=nfound+1;
        coords(nfound,:)=[ix iy]-n/2-1;
        amps(nfound)=amp;
    end;
    imacs(mximg2);
    title(nfound);
    drawnow;
    [amp ix iy]=max2di(mximg2);
end;
amps
return










for x=1:n
    for y=1:n
%         Theoretical values
        vecImage(:,:,x,y)=squeeze(vlist(:,ptrs(x,y),1:8));
        zr=vecImage(:,iGamma,x,y);
%         NCC
        zi=squeeze(ccs(x,y,:))/projamps(ptrs(x,y),iGamma);
        diff(x,y)=(zi-zr)'*(zi-zr); % residual power
        qi(x,y)=zi'*zi;  % image power
        qr(x,y)=zr'*zr;  % reference power (=1)
        qx(x,y)=zr'*zi;  % ncc
        zis(:,x,y)=zi;
        zrs(:,x,y)=zr;
    end;
end;
%%
subplot(221); imacs(qi);  title('Image power');
subplot(222); imacs(qx);  title('NCC');
subplot(223); imacs(m);
subplot(224); plot([sect(qi)  sect(qx) sect(diff)]);
subplot(224); plot([sect(qi)*1 sect(qx) sect(diff)]);

legend('i power','ncc','diff');

return

%%
% ------------------------------

zc=zeros(n,8);
y=129;
for x=1:n
    for g=4:4
        zc(x,g)=squeeze(vecImage(:,g,x,y))'*squeeze(ccs(x,y,:));
    end;
end;
    plot(zc);

disp('Assigned angles:');
ptr=zeros(nparts,1);
for i=1:nparts
    ptr(i)=ptrs(round(shifts(i,1)),round(shifts(i,2)));
    disp(round([angleList(ptr(i),:) (iGammas(i)-1)*dGamma ptr(i)]));
end;

x=round(shifts(1,1));
y=round(shifts(1,2));
disp('vectors:');
disp([vlist(:,ptr(1),iGamma) zrs(:,x,y) zis(:,x,y)]);
%     ptr

rv=reshape(ccrefs,n*n,nterms);
ishft=circshift(img,round(-shifts(1,:)));
imacs(ishft)
iv=reshape(ishft,n*n,1);
rv'*iv/(projamps(ptrs(x,y),iGamma))


    return

%% compare templates with reconstructions

subplot(221);
orig=allTemplates(:,:,ptr(i),iGamma);
imacs(orig);
subplot(222);
ei=reshape(eigenims,nt*nt,nterms);
ri=ei*vlist(:,ptr(i),iGamma)*projamps(ptr(i),iGamma);
recon=reshape(ri,nt,nt);
imacs(recon);
subplot(223);
plot([sect(orig) sect(recon)]);

oi=reshape(orig,nt*nt,1);
oexpansion=ei'*oi/projamps(ptr(i),iGamma);
disp([vlist(:,ptr(i),iGamma) oexpansion]);

%
orign=ifftshift(Crop(orig,n));
oin=reshape(orign,n*n,1);
cin=reshape(ccrefs,n*n,nterms);
oxn=cin'*oin/projamps(ptr(i),iGamma)




return
%%




load /Volumes/TetraData/Structures/AMPAR/3KG2RotMap5.8A.mat
nt=size(map,1);
map=map/64;  % approx amplitude correction (V-A scaling)
membraneOffset=-24/2;  % downsampled map by 2.

nAlpha=32; % about 10 degrees
nBeta=12;
nGamma=8;
symmetry=2;
nMirror=2;

[angleList inds]=rsListHemisphereAngles(nAlpha, nBeta);
nAngles=size(angleList,1);

% Add the gamma angles

dGamma=360/symmetry/nGamma;
allTemplates=zeros(nt,nt,nAngles,nGamma,nMirror);
angles=angleList;
angles(1,3)=0;  % add 3rd dimension
disp(['Making ' num2str(nGamma*nAngles*nMirror) ' templates']);
for j=1:2;  % upper, lower hemisphere
    if j==2
        mp=MirrorX(map);
    else
        mp=map;
    end;
    for i=1:nGamma
    gamma=(i-1)*dGamma;
    angles(:,3)=gamma;
    q=rsMakeTemplates(angles,mp);
    allTemplates(:,:,:,i,j)=q;
    end;
end;
%%
disp('Making eigenimages');
nterms=31;  % gets us to .9 of power.
ntotal=nterms;
        [eigenims vlist projamps termvar]=SphereMakeEigenRefs(allTemplates, nterms, ntotal);

        vlsPower=cumsum(vlist.^2);
        plot(reshape(vlsPower,nterms,nGamma*nAngles*nMirror));
        
        
%%        
n=256;
r0=100;

ptrs=rsGetAngleMap(n,r0,inds);
alphas=reshape(angleList(ptrs(:),1),n,n);
subplot(221);
imacs(alphas);
betas=reshape(angleList(ptrs(:),2),n,n);
subplot(222);
imacs(betas);
subplot(223);
imacs(ptrs);
drawnow;
%%
subplot(1,1,1);
 
for iGamma=1:16
img=zeros(n,n);
for i=1:32:256
    for j=1:32:256;
% % for ind=inds(1,:)
% %          for ind=[1 4 11 22 37 56 79 104 132 162 194 226];% 210 178 147 118 91 67 46 29 16 7 2 1];
        ind=ptrs(i,j);
%         disp([j ind angles(ind,:)]);
        q=allTemplates(:,:,ind,iGamma);
%         imacs(q);
%         title([i j]);
%         drawnow;
          tmp=ifftshift(Crop(allTemplates(:,:,ind,iGamma),n));
            img=img+circshift(tmp,[i j]);
    end;
end;
if iGamma==8 
    pause
end;
imacs(img);
title(iGamma);
drawnow;
end;
