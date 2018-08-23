% % RotatingDisplay.m
% 
% cd('/Users/fred/aEMCodeRepository/AMPAR')
% load KvMap.mat % loads 'map', 108 cubed; pixA=2;
% KvM=DownsampleGeneral(map,64,pixA/3.8);
% m0=rsRotateImage(KvM,40);
% m0=m0+.2*randn(64,64,64);
% figure(1);
% ImagicDisplay3(ExpandImage(m0(:,:,52:-1:8),4));

cd('/Users/fred/EMWork/Hideki/151117/KvLipo80slot3/Reconstructions/Recon96a2v/mrc')
[m2,s]=ReadMRC('i20av02.mrc');
% figure(2);
% ImagicDisplay3(ExpandImage(m2(:,:,53:-1:9),4));
% return

nGamma=20;
nBeta=10;
betaStep=6;
angs=zeros(nGamma*nBeta,3);
i=0;
for iBeta=1:nBeta
    beta=90+betaStep*(iBeta-1);
    for iGamma=1:nGamma
        gamma=iGamma*90/(nGamma+1);
        i=i+1;
        angs(i,:)=[0 beta gamma];
    end;
end;
ms=reMakeTemplates(m2, angs);

figure(6);
imovie(ms,.1);
