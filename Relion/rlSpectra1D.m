% rlSpectra1D.m
% From a star file, open the micrographs and compute the 1D power spectra.

% [nm,pa]=uigetfile('*.star');
% if isnumeric(pa)
%     return;
% end;
% 
% [nm,dat]=ReadStarFile([pa nm]);
d=dat{1};

nm=numel(d.rlnMicrographName);
sps=zeros(2048,nm);
for i=1:nm
    P=rlStarToCTFPars(d,i);
    [m,s]=ReadMRC(d.rlnMicrographName{i});
    disp(d.rlnMicrographName{i});
    sp2=fftshift(abs(fftn(m)).^2);
    sp0=Radial(sp2);
    pixA=1e4*d.rlnDetectorPixelSize(i)/d.rlnMagnification(i)
    sp1=RadialCTF(sp2,P,pixA);
%%
np=numel(sp1);
    asymp=sp1(np/2);
    sp1(1:20)=asymp;
    sp0(1:20)=asymp;
    np2=np/2;
    pts=1:np2;
    fs=(pts)/(2*np2*pixA);
    semilogy(fs,[sp0(pts) sp1(pts)]);
    title(i);
    drawnow;
    sps(:,i)=sp1;
end;
%%
    np2=np/2;
    pts=1:np2;
    fs=(pts)/(2*np2*pixA);

semilogy(fs,sps(pts,2:2:end));
%semilogy(fs,GaussFilt(sps(pts,2:2:end),.05,1));
xlabel('Spatial frequency, Å^{-1}');
ylabel('Spectral density');
semilogy(fs,sps(pts,1:2:end));
semilogy(fs,sps(pts,1:2:end));

save /ysm-gpfs/home/fjs2/scratch60/sps sps fs
    