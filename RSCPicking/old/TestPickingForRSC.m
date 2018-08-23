% TestPickingForRSC

load /Volumes/TetraData/Structures/AMPAR/3KG2RotMap5.8A.mat
nt=size(map,1);
map=map/64;  % approx amplitude correction (V-A scaling)
membraneOffset=-24/2;  % downsampled map by 2.

addNonParticles=1;

nAlpha=32; % about 10 degrees
nBeta=12;
nGamma=8;
symmetry=2;
nMirror=1;

[angleList inds]=rsListHemisphereAngles(nAlpha, nBeta);
nAngles=size(angleList,1);

% Add the gamma angles

dGamma=360/symmetry/nGamma;
allTemplates=zeros(nt,nt,nAngles,nGamma,nMirror);
angles=angleList;
angles(1,3)=0;  % add 3rd dimension
disp(['Making ' num2str(nGamma*nAngles*nMirror) ' templates']);
for j=1:nMirror;  % upper, lower hemisphere
    if j==2
        mp=MirrorX(map);
    else
        mp=map;
    end;
    for i=1:nGamma
    gamma=(i-1)*dGamma;
    tangles=angles;
    tangles(:,3)=gamma;
    q=rsMakeTemplates(tangles,mp);
    allTemplates(:,:,:,i,j)=q;
    end;
end;
size(allTemplates)
%%
if addNonParticles
% put in the non-particles
nrds=size(npRefs,4);
allTemplates(:,:,:,nGamma+1:nGamma+nrds,1)=npRefs*100;
size(allTemplates)
% templates are indexed as (x,y,ptrs,gamma,mirror)
else
nrds=0;
end;

%% CTF
c=abs(ifftshift(CTF(nt,3,.025,3,2,300,.07)));
allTemplatesCTF=reshape(allTemplates,nt,nt,nAngles*(nGamma*nMirror+nrds));
nim=size(allTemplatesCTF,3);
for i=1:nim
    m1=allTemplatesCTF(:,:,i);
    m2=real(ifftn(fftn(m1).*c));
    allTemplatesCTF(:,:,i)=m2;
end;

%%

%
disp('Making eigenimages');
nterms=21;  % gets us to .9 of power.
% nterms=12;  
ntotal=nterms;
        [eigenims vlist projamps termvar]=SphereMakeEigenRefs(allTemplatesCTF, nterms, ntotal);
        vlist=reshape(vlist,nterms,nAngles,nGamma*nMirror+nrds);
        vlist=vlist(:,:,1:nGamma);
        vlsPower=cumsum(vlist.^2);
        q=reshape(vlsPower,nterms,(nGamma)*nAngles*nMirror);
    subplot(223);
    plot(q);
    subplot(224);
    plot(projamps(:));
%%        
n=256;
r0=100;
figure(1); SetGrayscale;
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
figure(2);
SetGrayscale;
subplot(1,1,1);
 
for iGamma=1:8
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
% if iGamma==8 
%     pause
% end;
imacs(img);
title(iGamma);
drawnow;
end;


