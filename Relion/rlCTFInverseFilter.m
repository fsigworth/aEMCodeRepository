function m1=rlCTFInverseFilter(m,si,ds,compFraction,fHP)
% function m=rlCTFInverseFilter(m,si,ds,compFraction,fHP)
% Given the image stack m (downsampled by ds), boost the low frequencies below the 1st
% peak of the CTF.  If compFraction=1, the boost is the inverse of
% the ctf, up to a value of 10. fHP is the corner frequency of a Gaussian
% high-pass filter, in A^-1.

% if nargin<5
    fHP=0; % no high pass
% end;

k=.1; % limit on boost

[nx,ny,nim]=size(m);
n=[nx ny];
m1=m;
pixA=si.pixA*ds;
freqs=RadiusNorm(n)/pixA;  % frequencies of padded image
imi=0;
Hp=GaussHPKernel(n,fHP*pixA);

for i=1:nim
    if imi ~= si.miIndex(i)
        imi=si.miIndex(i);
        ct=Crop(si.ctfs(:,:,imi),n);  %%%???
        mi=si.mi{imi};
        [~,chi]=ContrastTransfer(freqs,mi.ctf(1));  % Just get chi
        peakMask=abs(chi)<0.5;  % mask up to first peak
        mxAmp=max2d(peakMask.*abs(ct));
        H=ones(n,'single');
        H(peakMask)=1-compFraction+(compFraction/mxAmp)./(abs(ct(peakMask))+k);
    end;
    m1(:,:,i)=real(ifftn(fftn(m(:,:,i)).*ifftshift(H.*Hp)));
end;

