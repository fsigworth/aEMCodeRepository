% TestImagePowerSpectrum
nExp=1;
ds=2;

filter='*.mat';
[fname fpath]=uigetfile(filter,'Select a micrograph info file');
cd(fpath);

load([fpath fname]);  % get the mi structure

disp('Base filename is');
disp(mi.baseFilename);

m=meReadImages(mi.rawFile);
%%
mp=mePreWhiten(m);
n=size(m,1);
freqs=(0:n/2-1)'/(n*mi.pixA);  % frequencies for evaluating the CTF
sp=zeros(n/2,3);
for i=1:3
    m1=meRemoveCrystal(mp(:,:,i),mi.pixA);
    sp(:,i)=RadialPowerSpectrum(m1)/mi.doses(i);
end;
%%

defs={};
for i=1:3
    defs{i}=num2str(mi.ctf(i).defocus);
end;

subplot(2,1,1);
loglog((freqs),sp);
xlabel('Spatial frequency');
ylabel('Power spectrum');
% axis([1e-3 inf .1 inf]);
legend(defs);
%%
cts=zeros(n/2,3);
for i=1:3;
    cts(:,i)=ContrastTransfer(freqs,mi.ctf(i));
end;
subplot(2,1,2);
semilogx(freqs,abs(cts));
% axis([0 inf 0 1]);
xlabel('Spatial frequency');
ylabel('abs(CTF)');
return

%%
if nExp==2
    [m effPixA]=ReadEMFile([mi.baseFilename 'em.mrc']);
elseif nExp==3
    [m effPixA]=ReadEMFile([mi.baseFilename 'ef.mrc']);
else % single exposure
    [m effPixA]=ReadEMFile([mi.baseFilename 'en.mrc']);
end;
m=RemoveOutliers(m);

m=Downsample(m,2048);  effPixA=effPixA*2;
    m=meRemoveCrystal(m,effPixA);
nd=size(m,1);

% effPixA=mi.pixA*ds;  % effective pixel size in our image
if nExp<3
    for i=1:numel(mi.ctf)
    mi.ctf(i).alpha=.05;
    end;
end;
freqs=Radius(nd)/(nd*effPixA);  % frequencies for evaluating the CTF
%%
% Compute the effective CTF after merging.  Note that this CTF is always
% positive; it doesn't include the contrast-reversal of the images.
if nExp>1
    [coeffs effctf]=meComputeMergeCoeffs(freqs, mi.ctf, mi.doses);
    flip=1;
else
    effctf=ContrastTransfer(freqs,mi.ctf(1));
    flip=-sign(effctf);
    effctf=abs(effctf);
end;
effctf=effctf.*CCDSqrtDQE(nd,ds);  % Correct for the camera's filtering

binf=2;

figure(1);
SetGrayscale;
subplot(2,2,1);
imacs(effctf);
subplot(2,2,2);
plot(sectr(freqs),sectr(effctf));
xlabel('spatial frequency, A^{-1}');
subplot(2,2,3);
imacs(BinImage(m,binf));

%%
% Find the first and second zeros
thresh=.01;
radialctf=sectr(effctf);
pt=find(radialctf<thresh);
pt1=pt(1);
st=find(diff(pt)>1);
pt2=pt(st(2));
% find the maximum between the first and second zeros.
% We'll normalize the filtered image to this level.
mx2=max(radialctf(pt1:pt2));
mx2m=mx2/2;  % LF amplitude
R=Radius(nd);
HLow=ones(nd,nd);
ptLow=(effctf>mx2m | R<pt1/2) & (R< pt1);
HLow(ptLow)=mx2m./(effctf(ptLow));
ectf2=HLow.*effctf;
subplot(2,2,2);
plot(sectr(freqs),[sectr(effctf) sectr(ectf2)]);
ectfs(:,nExp)=sectr(effctf);
ectfs2(:,nExp)=sectr(ectf2);

HLow=HLow.*flip;
mf=real(ifftn(fftn(m).*fftshift(HLow)));
subplot(2,2,4);
imacs(BinImage(mf,binf));

figure(2);
clf;
SetGrayscale;
imacs(BinImage(mf,1));

imwrite(uint8(imscale(mf)),[mi.baseFilename 'mcw' num2str(nExp) '.tif']);

%%
figure(3);
defstr='defocus = ';
for i=1:3
subplot(3,1,i);
plot(sectr(freqs),[ectfs(:,i) ectfs2(:,i)]);
ylabel('CTF');
defstr=[defstr '    ' num2str(mi.ctf(i).defocus)];
title(defstr);
end;
xlabel('Spatial frequency, A^{-1}');
