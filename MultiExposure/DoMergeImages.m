% function [mc effctz doses Tmats ctfs]= meMergeImages(fnames,defoci)
% function [minfo]= meMergeImages(minfo,defoci)
% Merge a defocus series (2 or more images)

% This presently runs as a script; eventually will be a function called in
% a loop with different sets of filenames.
% Fixed-pattern compensation is now done.
% The "output" variables are the merged image mc, the effective CTF effctz,
% individual image contributions mcs, the doses, and the affine
% transformation matrices Tmats.

% The major function calls used here are
%
% meReadImages: reads, removes outliers, windows and computes doses
% mePreWhiten: Operates on all images with the CCD pre-whitening filter
% meSetCTFitPars: Initializes the CTF-fitting parameters
% ctfit2: main CTF fitting function
% meAlignMultiExposures: Alignment of the images
% meCombineImages: compute the optimum weightings and merge the images.

tic

disp(' ');
disp('meMergeImages...');

% pick up the path
pa='/Users/fred/EMWork/Liguo/10sep19a/';
% pa='/Volumes/fred/EMWork/Liguo/10sep19a/';
% pa='/Users/fred/EMWork/Liguo/10dec18a/';
d=dir(pa);
fnames=[];

nfiles=3;
fileoffset=1;
xoffset=2;
for i=1:nfiles
    fnames{i}=[pa d(2+xoffset+2*i+6*fileoffset+2*(3-nfiles)).name];
end;
char(fnames)
[pa1 basename imtype ext]=ParseLeginonName(fnames{1});
basename
%%

% concatenate file names
% fnames={[pa '10sep19a_a_00006gr_00017sq_v01_00021hl_v02_00004em.mrc'],...
%         [pa '10sep19a_a_00006gr_00017sq_v01_00021hl_v02_00004en.mrc']};
%     ,...
%         [pa '10sep19a_a_00006gr_00017sq_v01_00021hl_v02_00004en.mrc']};
disp(char(fnames));

y=0;
m=[];
DoRemoveSpots=1;  % 1 for removing crystal-lattice periodicity before running ctfit

if nfiles==3
defoci=[10 3 1];  % defocus values in um
else
    defoci=[3 1];  % defocus values in um
end;

disp('meReadImages');
[m pixA doses]=meReadImages(fnames);
nim=size(m,3);

% operate with the pre-whitening filter
disp('mePreWhiten');
m=mePreWhiten(m);

%%
% Make sure the defocus values are sorted from lowest to highest;
% if not, re-order defoci, the m stack and the file names.
if ~issorted(defoci)
    %%
    [defvals inds]=sort(defoci);
    %%
    mtemp=m;
    dtemp=doses;
    tnames=fnames;
    for i=1:size(m,3)
        mtemp(:,:,i)=m(:,:,inds(i));
        tnames{i}=fnames{inds(i)};
        doses(i)=dtemp(inds(i));
    end;
    m=mtemp;
    mtemp=[];  % deallocate this big array
    defoci=defvals;
    fnames=tnames;
end;
%%
% Create starting CTF fitting values
CTFitPars=meSetCTFitPars(defoci,pixA);

% Do the CTF fitting
n=size(m,1);
for i=1:nim
    mc=m(:,:,i);
    if defoci(i)>=1.5 % large defocus; fit CTF with binned image
        ds=2;
        maxres=10;
        minres=100;
    else            % small defocus: fit CTF with full image
        ds=1;
        maxres=5;
        minres=50;
    end;
    ns=n/ds;
    [mc fmc]=Downsample(mc,ns);  % reduce image size by ds
    
    figure(4); clf;
    SetGrayscale;
    imacs(mc);
    drawnow;
    if DoRemoveSpots
        fmcs=fftshift(fmc);
        ds1=max(ds,2);
        nr=n/ds1;
        if ds1>ds
            fmcs=Crop(fmcs,nr);
        end;
        minr=n*pixA/70;
        ThreshSD=5;
        display=1;
        disp('RemoveSpots');
        [spc pmask]=RemoveSpots(abs(fmcs).^2,minr,ThreshSD,display);
        mc=real(ifftn(fftshift(1-Crop(1-pmask,ns)).*fmc));
    end;
    
    disp(['CTF fitting ' num2str(i)]);
    [P c]=ctfit2(mc,CTFitPars(i),pixA*ds,maxres,minres);
    
    P.res=pixA;  % Change it back to the original pixel size
    CTFitPars(i)=P;
end;
%% reclaim space
pmask=[];
spc=[];
fmc=[];
fmcs=[];
mc=[];
%%
AlignDS=4;

% Align the images
Tmats=meAlignMultiExposures(m,pixA,CTFitPars,doses,AlignDS);

%%
% Combine the images
[mc effctz mcs]=meCombineImages(m,pixA,CTFitPars,Tmats,doses,ds);
mc=-mc;  % go back to reversed contrast

%%
% Wiener filter
fmc=fftshift(fftn(mc));
k=.01;
mcw=real(ifftn(fftshift(fmc.*effctz./(k+effctz.^2))));
fc=.08;
% Additional sharp filter
% mcw=GaussFilt(fmcw,fc*pixA);
mcw=SharpFilt(mc,fc*pixA,fc*pixA/4);
imacs(mcw);
imwrite(uint8(imscale(mcw)),['x' basename 'mf2.tif'],'tiff');

% return

%%
%  Remove crystal pattern from merged image
ns=n/ds;
sp=Crop(fftshift(abs(fftn(mcw)).^2),n/4);
[spc pmask]=RemoveSpots(sp,.8*4096*pixA/58,5,1);
%
msk=(fuzzymask(n/4,2,300,10));
pmaskc=msk.*pmask+(1-msk);
imacs(pmaskc); drawnow;
fmask=1-Crop(1-pmaskc,ns);  % fourier mask for image.
fm1=fftshift(fftn(mcw));
mcm=real(ifftn(fftshift(fm1.*fmask)));
subplot(1,1,1); imacs(mcm);
imwrite(uint8(imscale(mcm)),['x' basename 'mf2s.tif'],'tiff');
pixA*ds

toc
