% rsTestFourierReconstruction
% Test the rsFourierReconstruction function

n=64;
disp('Loading model');
load 3KG2RotMap2.9A.mat
n0=size(map,1);
m0=Crop(map,128);
m1=Downsample(m0,n);
figure(2);
ShowSections(m1);
lambda=EWavelength(200);
sigmaN=500;  % noise SD
k=.1;  % Wiener constant

% ri.pixA=3*128/n;
% ri.n             =  n;  % box size
% ri.nSteps= [36 12 12];  % Alpha, beta, gamma steps
% ri.symmetry = 2;  % We're assuming C symmetry only
% ri.gammaMode='variable';
% 
% angles=rsMakeTemplateAngles3(ri);

ri.pixA=3*128/n;
ri.n             =  n;  % box size
ri.nSteps= [0 12 12];  % Alpha, beta, gamma steps
ri.nSteps= [0 24 24];  % Alpha, beta, gamma steps
ri.symmetry = 2;  % We're assuming C2 symmetry only
ri.gammaMode='variable';

angles=rsMakeTemplateAngles3(ri);





disp('Making templates');
nangs=size(angles,1)
templates=rsMakeTemplates(angles,m1);

%%
fnorms=(zeros(n,n,nangs));
fprojs=(zeros(n,n,nangs));
imgs=(zeros(n,n,nangs));
disp('Filtering');
for i=1:nangs
    def=rand*2+1;  % defocus between 1 and 3
    B=50+50*def;
    c=abs(CTF(n,ri.pixA,lambda,def,2,B,.2));  % Assume phase-flipping
%     c=-(CTF(n,ri.pixA,lambda,def,2,B,.2));  % no phase-flipping
    noise=sigmaN*randn(n,n);
    img=real(ifftn(fftn(templates(:,:,i)+noise).*ifftshift(c)));
    fprojs(:,:,i)=c.*fftshift(fftn(ifftshift(img)));  % filter with a single ctf
    fnorms(:,:,i)=c.^2;
    imgs(:,:,i)=img;
end;
disp('Reconstructing');
[vol rv0]=rsFourierReconstruction(angles, fprojs, fnorms, ri.symmetry, k);
disp('Done.');
%
% Show coverage of the Fourier volume
figure(2);
ShowSections(sqrt(abs(fftshift(fftn(rv0)))));
figure(1);
ShowSections(vol.*fuzzymask(n,3,n*0.45,n*.05));


