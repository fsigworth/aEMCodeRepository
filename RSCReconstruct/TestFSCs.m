% testFSCs.m

cd /Users/fred/EMWork/Hideki/140625/KvBetaLiposome_new_pH8_KvLipo11_slot1/Reconstructions/Recon128
load iterVols.mat

%%
niters=20;
fs=zeros(47,niters);
for ind=1:20;
    v1=outputVols1(:,:,:,ind);
    v2=outputVols2(:,:,:,ind);
    msk=fuzzymask(96,3,40,10);
%     msk=fuzzymask(96,3,20,10,[49 49 38]);
    % msk=1;
    figure(5);
    ShowSections2(v1.*msk);
    figure(6);
    ShowSections2(v2.*msk);
    subplot(3,3,9);
    fs(:,ind)=FSCorr(v1.*msk,v2.*msk);
    plot(fs);
    drawnow;
end;
%%

% Rotationally align
[x,y,z]=ndgrid(-48:47);
r=sqrt(x.^2+y.^2+z.^2);
% rs=ifftshift(r.^2);
% cc3=fftshift(real(ifftn(fftn(v1.*msk).*rs.*conj(fftn(v2.*msk)))));
% figure(6);
% ShowSections2(cc3);
%
% [val p]=max3di(cc3)


%%
ng=11;
cr0=v1(:)'*v2(:);
cr=zeros(ng,1);
fs=zeros(47,ng);
for i=1:ng
    gamma=-5.5+.5*i;
    v1r=rsRotateImage(v1,gamma*ones(96,1));
    cr(i)=v1r(:)'*v2(:);
    disp([gamma cr(i)])
    fs(:,i)=FSCorr(v1r,v2);
    plot(fs);
    drawnow;
end;
