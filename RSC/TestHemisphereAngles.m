% TestHemisphereAngles.m

nAlpha=64;
nBeta=32;
n=127;
r0=60;

[angleList inds]=rsListHemisphereAngles(nAlpha, nBeta);
ptrs=rsGetAngleMap(n,r0,inds);
imacs(ptrs)
alphas=reshape(angleList(ptrs(:),1),n,n);
subplot(221);
imacs(alphas);
betas=reshape(angleList(ptrs(:),2),n,n);
subplot(222);
imacs(betas);