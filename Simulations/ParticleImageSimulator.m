% ParticleImageSimulator.m

% The angles are phi (spin about the z axis), theta (tilt) and psi
% (rotation in image plane). So here I spin the particle, then lay it on
% its side, then spin in plane.

% It takes a long time to project volumes larger in size than say 128.
% My code requires that the size be a multiple of 8.

n=128; % volume/image size

[map,s]=ReadMRC('emd_21865Ds2.mrc');

m=DownsampleGeneral(map,n);
pixA=s.pixA*size(map,1)/n;
disp(['Pixel size is ' num2str(pixA) 'A']);
% This is the pixel size, about 2A

angles=[ 0  0  0
    10 0 0
    20 0 0
    30 0 0
    30 45 0
    30 90 0
    30 90 10
    30 90 20
    30 90 30];

projs=rlMakeTemplates(angles,m);

% CTF parameters
lambda=EWavelength(300);
def=2; % 2um defocus
Cs=2.7; % typical for Krios
B=60; % very good micrographs
alpha=.07; % amplitude contrast

c=CTF(n,pixA,lambda,def,Cs,B,alpha); % no astigmatism parameters given
ci=ifftshift(c); % zero frequency shifted to origin
nim=size(projs,3);
cProjs=zeros(n,n,nim,'single');
for i=1:nim
    cProjs(:,:,i)=real(ifftn(fftn(projs(:,:,i)).*ci));
end;

subplot(2,2,1);
imags(projs(:,:,1)); % Show the first one
title('First projection');

subplot(2,2,2);
imags(c); % Show the CTF, zero frequency in the center
title('CTF');

subplot(2,2,3);
imags(cProjs(:,:,1)); % CTF-filtered first projection
title('First CTF-filtered');

sigma=2;
sProjs=cProjs+sigma*randn(n,n,nim);
subplot(2,2,4);
imovie(sProjs,.2); % 'movie' of all.




