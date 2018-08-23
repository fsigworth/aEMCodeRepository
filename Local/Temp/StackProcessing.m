% How I'd process the images, yielding the final array filtImgs
% -Fred
stk=ReadMRC(fileName);
% alternatively, with startSlice starting at 1:
% stk=ReadMRC(fileName,startSlice,numSlices);

n=size(stk,1);  % Had better be even, or my functions will have problems
nim=size(stk,3);

pixA=6.35/60765.550200*10000;  % pixel size in angstroms
lambda=EWavelength(300);
ampContrast=.1;  % Their given value
B=1;  % not sure what to use for B-factor, maybe 30 in units of A^2

% Estimate noise and signal spectra
msk=fuzzymask(n,2,0.45*n,.05*n);  % soft mask to separate rim from center
% The following call to RimSpectrum will take a long time to compute if the stack is large.
noiseSpec1d=RimSpectrum(stk,1-msk);  % get the average 1d spectrum using the rim of the imageswide on the sides)

noiseSpec2d=ifftshift(ToRect(noiseSpec1d,n));  % Polar to rect: in this case, makes it circularly symmetric
sigSpec1d=RadialPowerSpectrum(stk,1);
sigSpec1d(round(0.4*n):end)=mean(SigSpec1d(round(0.3*n):round(0.4*n)));  % force it to be constant up to Nyquist
sigSpec2d=ifftshift(ToRect(sigSpec1d,1));

k=noiseSpec2d./sigSpec2d;  % Or whatever, your Wiener 'constant'
   
% Do Wiener filtering
% I'm assuming that you already have initialized d=vector of defocus values
% in microns (the given values in the file look like they are in angstroms)
filtImgs=zeros(size(stk),'single');
for i=1:nim
    c=ifftshift(CTF(n,pixA,lambda,d(i),0,B,ampContrast));  % This gives zero frequency at (1,1)
    img=stk(:,:,i);
    filtImgs(:,:,i)=real(ifftn(fftn(img).*c./(k+c.^2)));
end;
