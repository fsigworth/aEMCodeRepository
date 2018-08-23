function [msub, minfo]=meSubtractVesicles(m, minfo)
% function [minfo msub]=meSubtractVesicles(m, minfo);
% Given an image m and the info structure minfo, find and subtract vesicles
% and return the updated minfo containing the vesicle coordinates, and the
% subtracted image.

% serverPath='/Volumes/TetraData/EMWork/Liguo/FavoriteBK/';
% infoFilename=[serverPath '10dec18aMi/10dec18a_bkoWA15B_00054gr_00043sq_v01_00007hl_v01_00004minfo.mat'];
% load(infoFilename);
% mi.serverPath=serverPath;
% minfo=mi;
% 
% m=ReadEMFile([mi.serverPath mi.localPath mi.baseFilename 'merge.mrc']);

displayOn=1;
CCDCountsPerE=14;


n=size(m,1);
ds0=minfo.nPix(1)/n;  % downsampling of the image we're working with.
pixA0=minfo.pixA*ds0;  % pixA of merged image
targetPixA=10;  % maximum pixel size
ns=NextNiceNumber(n*pixA0/targetPixA,5,4);  % multiple of 4, max factor 5.
if ns<n
    ms=Downsample(m,ns);
end;
fms=fftn(ms);
ds=n/ns;  % downsampling factor
pixA=ds*pixA0;  %pixA in the downsampled image.

mbnThickness=50/pixA;  % membrane thickness in pixels for model
rmin=100/pixA; % minimum liposome radius in pixels
rmax=400/pixA; % maximum liposome radius
rstep=10/pixA; % radius step
nbin=2;  % binning of image for display
mbnPhase2=2*1.2e-3;  % radians per A x 2 for phase contrast
voltPhase2=.73e-3; % radians per V*A x 2 for phase contrast


figure(1);
SetGrayscale;
subplot(2,3,1);
imacs(ms);

% Get the effective CTF from the merging.
freqs=Radius(ns)/(ns*pixA);
[~, effctz]=meComputeMergeCoeffs(freqs, minfo.ctf, minfo.doses);
H=fftshift(effctz.*CCDSqrtDQE(ns,ds));  % CTF filter

%%
disp('Constructing references...');

nrsteps=round((rmax-rmin)/rstep+1);
%
rrefs=single(zeros(ns,ns,nrsteps));
frefs=single(zeros(ns,ns,nrsteps));
rads=zeros(nrsteps,1);

subplot(2,3,2);
nd=NextNiceNumber(3*rmax/pixA);  % size of box to display
nd=min(nd,ns);

for i=1:nrsteps
    r=rmin+(i-1)*rstep;  % radius in pixels
    rads(i)=r;  % radius in A
    v = VesicleDensity(ns,r,mbnThickness);
    fv=-H.*fftn(fftshift(v));
    frefs(:,:,i)=fv;
    rv=fftshift(real(ifftn(fv)));
    refPower(i)=rv(:)'*rv(:);
    rrefs(:,:,i)=rv;
    if displayOn
        imacs(rv);
        title(i);
        drawnow;
    end;
end;

% 
% %     v*pixA*mbnPhase2;  % Theoretical contrast for a vesicle
%     
%     imacs(Crop(rv,nd));
%     title(['Radius ' num2str(r*pixA)]);
%     drawnow;
%     
%     fs1=fftn(fftshift(rs(:,:,i)));  % FT of sphere
%     fv1=fftn(fftshift(rv(:,:,i)));        % FT of vesicle
%     fv(:,:,i)=fv1;  % Fourier-transformed reference for cross-correlations
%     fs(:,:,i)=fs1;
%     fv1f=fv1.*F;    % High-pass filtered FT reference
%     fs1f=fs1.*F;
%     sv(i)=abs(fv1(:)'*fv1f(:))/n^2; % normalization of amplitude.
%     ss(i)=abs(fs1(:)'*fs1f(:))/n^2;
% end;

%%
% Compute all the cross-correlations, using a weight of r^-2 to equalize
% the CC peak values
disp('Computing cross-correlations');
ccs=single(zeros(ns,ns,nrsteps));
for i=1:nrsteps
    ccs(:,:,i)=rads(i).^-2*real(ifftn(fms.*conj(frefs(:,:,i))));
end;

[ccmx ccmi]=max(ccs,[],3);
subplot(2,3,2);
imacs(ccmx);
drawnow;

%%
disp('Finding vesicles');


ccmx2=ccmx;
ms2=ms;

% Find the first peak.
[presentmax i j]=max2di(ccmx2);
[interpmax rind]=max1di(ccs(round(i),round(j),:));

% scan all the CC peaks
maxfract=0.4;  % Minimum CC value accepted, relative to the global max

coords=zeros(2,1);
amps=zeros(2,1);
radii=zeros(2,1);
refis=zeros(2,1);
model=zeros(ns,ns);
umodel=zeros(ns,ns);
globalmax=presentmax;
nfound=0;

while presentmax>globalmax*maxfract
    nfound=nfound+1;
    % find out which reference gave rise to the maximum
    coords(:,nfound)=[i;j]-1;  % origin is (0,0)
    [mxv refii]=max1di(ccs(round(i),round(j),:));
    
    refi=ccmi(round(i),round(j));
%     [refi refii]
    refis(nfound)=refi;        % Which reference matched the best
%     refr=rads(refi);            % liposome radius in pixels
    refr=rmin+(refii-1)*rstep;
    radii(nfound)=refr;
    ampi=presentmax*rads(refi)^2/refPower(refi); % undo the we ighting by r^-2 in the CC
    amps(nfound)=ampi;
    
    %     ampi=12.4; %  fix the median value.
    % simple integer-shifted reference
    %      sref=ampi*circshift(rv(:,:,refi),round([i j])-n/2-1);
    %      sref=real(ifftn(fftn(sref).*PF));
    vref=ampi*circshift(rrefs(:,:,refi),round([i j]-1));
    umodel=umodel+fftshift(vref);  % unfiltered model
    msub=ms-umodel;
    % Blank the region in the cross-correlation
    blank=fuzzymask(ns,2,refr+mbnThickness/2,mbnThickness,[i j]);
    ccmx2=ccmx2.*(1-blank);
%     subplot(2,3,4); imacs(ms2);
%     title(nfound);
%     subplot(2,3,5); imacs(ccmx2);
    %      subplot(2,3,4); imacs(BinImage(model,nbin));
%     drawnow;
    [presentmax i j]=max2di(ccmx2);
    
end;
    subplot(2,3,4); imacs(msub);
    title(nfound);
    subplot(2,3,5); imacs(ccmx2);
    subplot(2,3,2); imacs(umodel);
    subplot(2,3,6); imacs(ccmx);
subplot(2,3,3);
plot(radii,amps,'k.');
xlabel('Vesicle radius, Å');
ylabel('Image amplitude');
drawnow;
% %%
% % Create new vesicle model
% model=zeros(ns,ns);
% for i=1:nfound
%     model=model+amps(i)*VesicleDensity(ns,radii(i),mbnThickness,coords(:,i)+1);
% end;
% subplot(2,3,2); imacs(-model);

%% Insert vesicle coordinates into info file
minfo.vesicle.x=single(zeros(nfound,1));
minfo.vesicle.y=single(zeros(nfound,1));
minfo.vesicle.r=single(zeros(nfound,1));
minfo.vesicle.amp=single(zeros(nfound,1));
for i=1:nfound
    minfo.vesicle.x(i)=coords(1,i)*ds0*ds;  % in pixels of original image
    minfo.vesicle.y(i)=coords(2,i)*ds0*ds;
    minfo.vesicle.r(i)=radii(i)*ds0*ds;
    minfo.vesicle.amp(i)=amps(i);
end;

minfo.vesicleModel=mbnThickness*ds0*ds;  % the single scalar is unform mbn thickness

% mi=minfo;
% save(infoFilename,'mi');
