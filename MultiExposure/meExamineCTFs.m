% meExamineCTFs
% ds=2;
nExp=3;

filter='*.mat';
[fname fpath]=uigetfile(filter,'Select a micrograph info file');
cd(fpath);

load([fpath fname]);  % get the mi structure

disp('Base filename is');
disp(mi.baseFilename);
raw=0;
fname=[mi.procPath mi.baseFilename];
    if FileExists([fname 'mc.mrc'])
        fname=[fname 'mc.mrc'];
    elseif FileExists([fname 'm.mrc'])
        fname=[fname 'm.mrc'];
    elseif FileExists([mi.imagePath mi.baseFilename 'en.mrc'])
        fname=[mi.imagePath mi.baseFilename 'en.mrc'];
        raw=1;
    else
        error('Can''t find an image file');
    end;
    [m effPixA]=ReadEMFile(fname);
    if raw
        m=RemoveOutliers(m);
        m=Downsample(m,2048);  effPixA=effPixA*2;
    end;
nd=size(m,1);

% effPixA=mi.pixA*ds;  % effective pixel size in our image
% if nExp<3
%     for i=1:numel(mi.ctf)
%     mi.ctf(i).alpha=.05;
%     end;
% end;
freqs=Radius(nd)/(nd*effPixA);  % frequencies for evaluating the CTF
%%
% Compute the effective CTF after merging.  Note that this CTF is always
% positive; it doesn't include the contrast-reversal of the images.
if ~raw
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
thresh=.001;
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
imacs(GaussFilt(mf,.1));

imwrite(uint8(imscale(mf)),[mi.baseFilename 'mcw' num2str(nExp) '.tif']);

%%
figure(3);
defstr='defocus = ';
for i=1:size(ectfs,2)
subplot(3,1,i);
plot(sectr(freqs),[ectfs(:,i) ectfs2(:,i)]);
ylabel('CTF');
defstr=[defstr '    ' num2str(mi.ctf(i).defocus)];
title(defstr);
end;
xlabel('Spatial frequency, A^{-1}');
