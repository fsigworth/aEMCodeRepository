% TestABReconstruct

%% Load Example Inputs from Fred's Data Structure and Convert to Hemant's Notation
pa=ParsePath(which('adaptiveBasisReconstructionFRED'));
pa=[pa 'dataFromFred/'];
dat=load([pa 'AlpsGoodies.mat']);
sigmaN=30;

[n, ny, nim]=size(dat.imgsOrig);
ctf0=abs(CTF(n,5,.025,1,2,10,.07));
% ctf0=ones(n,n,'single');
amps=repmat(rand(1,1,nim),n,n);
ctfs=repmat(ctf0,1,1,nim).*amps;
ctfs0=repmat(ifftshift(ctf0),1,1,nim);  % zero-origin
imgs=amps.*(real(ifft2(fft2(dat.imgsOrig).*ctfs0))+sigmaN*randn(n,n,nim));
vol=abReconstruct(imgs,ctfs,dat.allSimAngles);
figure(5);
ShowSections2(vol);
fsc=FSCorr(vol,dat.origVols);
figure(4); clf
plot(fsc);
figure(2);
ImagicDisplay2(imgs);
