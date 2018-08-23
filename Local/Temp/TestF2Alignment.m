% TestF2Alignment.m
% Examine the fixed-pattern noise in Falcon II movies.

figure(1);
SetGrayscale;


cd('/Volumes/WD2Blue/EMWork/Hideki/140722/KvBetaLiposome_sol_slot1/movie_frames/sq05');
ind=4;
d=dir;
name=d(ind).name
[m,pixA]=ReadEMFile(name);
mf=single(m);
ms=sum(mf,3);
imacs(ms);
imgMean=mean(ms(:))
%%
cc17=fftshift(real(ifftn(fftn(mf(:,:,1)).*conj(fftn(mf(:,:,6))))));
ccSm=Crop(cc17,64);
ccCtr=Crop(cc17,7);
q=(ccSm(31:35,31:35)-ccSm(31,31))/-2e8
figure(1);
plot(ccSm(26:38,31:35),'.-');
% imacs(del2(ccSm));
%% Try removing the artifact by forcing the Laplacian to be continuous
for i=1:3  % iterate to remove errors
    d2=del2(ccCtr);
    
    % Define the points from which we compute the average Laplacian
    msk=[0 0 0 0 0 0 0
        0 0 0 1 0 0 0
        0 0 1 1 1 0 0
        0 0 1 0 1 0 0
        0 0 1 1 1 0 0
        0 0 0 1 0 0 0
        0 0 0 0 0 0 0];
    dm=msk(:)'*d2(:)/(msk(:)'*msk(:));  % compute the mean of 2nd derivative
    dc=d2(4,4);
    ccCtr(4,4)=ccCtr(4,4)+1*(dc-dm);
end;
ccSm1=ccSm;
ccSm1(33,33)=ccCtr(4,4);
% plot([sect(ccSm) sect(ccSm1)]);
figure(2);
plot(ccSm1(26:38,31:35),'.-');
% imacs(del2(ccSm1));

q1=(ccSm1(31:35,31:35)-ccSm1(31,31))/-2e8
ccSm=ccSm1;
