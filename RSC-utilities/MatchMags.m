% MatchMags.m
% Try matching magnifications between K2 and F2 data

% Get the K2 data
cd('/Volumes/fred/EMWork/Hideki/140625/KvBetaLiposome_new_pH8_KvLipo11_slot1/Reconstructions/Recon96j/mrc')
[v1,s1]=ReadMRC('i19av01.mrc');
v1=v1+ReadMRC('i19bv01.mrc');
figure(1);
ShowSections2(v1,[],45);
title('v1');

% Get the F2 data
cd('/Volumes/IndyT3/data/150119/KvLipo14slot3/Reconstructions/Recon96o/mrc')
[v9,s9]=ReadMRC('i19cv01.mrc');
figure(2);
ShowSections2(v9,[],45);
title('v9');

n=size(v1,1);

msk=fuzzymask(96,3,24,8,[49 49 39]);

% msk=1

startingRelMag=s9.pixA/s1.pixA  % F2 mag is smaller

nm=9;
dm=.02;
m0=1;
mags=m0+dm*(-nm/2:nm/2-1);
n3=n*[1 1 1];
mag=Simplex('init',m0,dm);
niters=25;
mags=zeros(niters,1);
ccfs=zeros(niters,1);
for iter=1:niters
    v9m=DownsampleGeneral(v9,n,mag);
    v9m=circshift(v9m,[0 0 5]);
    % Align the volumes
    cct=fftshift(real(ifftn(fftn(v9m).*conj(fftn(v1)))));
    [val,coords]=max3di(cct);
    coords=coords-ceil((n+1)/2);
    P=FourierShift(n3,-coords);
    v9ms=real(ifftn(fftn(v9m).*P));
    
    figure(2);
    ShowSections2(v9ms,[],45);
    title([iter mag]);
    drawnow;
    v9msk=v9ms.*msk;
    denom=sqrt((v1(:)'*v1(:))*(v9msk(:)'*v9msk(:)));
    ccf=v1(:)'*v9msk(:)/denom;
    disp([mag ccf]);
    mags(iter)=mag;
    ccfs(iter)=ccf;
    mag=Simplex(-ccf);
end;
%%
figure(3);
plot(mags,ccfs,'ko');

fsc=FSCorr(v9m,v1);
df=1/(n*s1.pixA);
fs=df:df:numel(fsc)*df;
figure(4);
plot(fs,[fsc 0*fsc]);
