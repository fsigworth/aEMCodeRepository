% rsTestGriddingReconstruction
% Test the rsFourierReconstruction function

n=64;
n1=64;

disp('Loading model');
load 3KG2RotMap2.9A.mat
n0=size(map,1);
m0=Crop(map,128);
m1=Downsample(m0,n);
figure(2);
ShowSections(m1);
drawnow;
lambda=EWavelength(200);
sigmaN=500;  % noise SD
k=.001;  % Wiener constant

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
norms=(zeros(n,n,nangs));
imgs=(zeros(n,n,nangs));
disp('Filtering');
for i=1:nangs
    def=rand*2+1;  % defocus between 1 and 3
    B=50+50*def;
    c=abs(CTF(n,ri.pixA,lambda,def,2,B,.2));  % Assume phase-flipping
    noise=sigmaN*randn(n,n);
    img=real(ifftn(ifftshift(c).*fftn(templates(:,:,i)+noise)));
    norms(:,:,i)=fftshift(real(ifftn(ifftshift(c.^2))));
    imgs(:,:,i)=real(ifftn(ifftshift(c).*fftn(img)));  % operate again by ctf
end;

%%
k=.01;  % Wiener constant

disp('Reconstruction');
[vol norm volf]=rsGriddingReconstruction(angles,imgs,norms,k);
cvol=Crop(vol,n1);
cnorm=Crop(norm,n1);  % crop to smaller size for display
figure(1);
fnorm=fftshift(abs(fftn(norm)));
%
ShowSections(sqrt(fnorm));
figure(2);
ShowSections(cvol);
subplot(3,3,8);
imagesc(img);
axis equal off;
xlabel('input projection');
