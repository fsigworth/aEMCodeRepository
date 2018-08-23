% ComputeTomo.m
% load and align a tomogram stack.

ds1=2;
ds2=2;
ds=ds1*ds2;

doAutoAlign=0;
doManualAlign=0;

% ----------------Get the stack files---------------

[fname, pa]=uigetfile('*.st.mdoc','Select a mosaic doc file');
if isnumeric(pa)  % user clicked Cancel
    return
end;

cd(pa);

%%


disp(['reading ' fname]);
%%
mdoc=fopen(fname);
lines={};
i=0;
j=0;
zVals=[];
assigned=false(0,2);
s=fgetl(mdoc);

% Get image info
mosaic=zeros(0,1,'single');

while ~feof(mdoc) && numel(s)>0
    c=textscan(s,'%s = %f %f','emptyvalue',0);
    tok=char(c{1});
    switch tok
        case 'PixelSpacing'
            pixelSpacing=c{2};
    end;
    s=fgetl(mdoc);
end;
%  a blank line follows the header.
nim=0;
angles=zeros(0,1);
while ~feof(mdoc)
    if numel(s)>0
        c=textscan(s,'%s = %f %f %f','emptyvalue',0);
        tok=char(c{1});
        switch tok
            case 'TiltAngle'
                nim=nim+1;
                angles(nim)=c{2};
        end;
    end;
    s=fgetl(mdoc);
end;

angles'

cosines=cosd(angles)';

stackName=fname(1:end-5);
disp(['Reading ' stackName]);
[st, pixA]=ReadEMFile(stackName);
fst=st;
[nx, ny, nim]=size(st);
nxw=floor(nx/ds);
nyw=floor(ny/ds);
n0=NextNiceNumber(max(nx,ny));

%%
disp('Removing outlers');
thresh=Percentile(fst,1-1e-7);
fst(fst>thresh)=thresh;

%%
disp('Binning');
stbin=zeros([n0/ds1 n0/ds1 nim],'single');
for i=1:nim
    q=single(st(:,:,i));
    me=mean(q(:));
    stbin(:,:,i)=SharpFilt(BinImage(Crop(q,n0,0,me),ds1),.5/ds2,.1/ds2);
end;

%  Auto-alignment
% Find the zero-tilt image
iz=find(abs(angles)<1,1,'first');

%% Attempt to align the images
if doAutoAlign
    nw=size(stbin,1);
    kp=2;
    pars=zeros(nim,4);
    ctr=nw/2+1;
    ix=ctr;
    iy=ctr;
    niters=1;
    fcc=.2;
    shStack=stbin;
    % for i1=iz-1:-1:1
    ist=iz+1;
    ien=nim;
    ivals=[iz-1:-1:1 iz+1:nim];
    for i1=ivals
        if i1==iz+1 % Shift the center
            ix=ctr+300;
            iy=ctr-200;
        end;
        prior=kp*RadiusNorm(nw,[ix iy]).^2;
        p0=[-5.2 1];
        p=Simplex('init',p0,[.1 .01],[0 1]);
        p=p0;
        for i=1:niters
            theta=p(1);
            ang1=angles(i1)*p(2);
            T0=[cosd(ang1)/cosines(iz) 0 0; 0 1 0; 0 0 1];
            R0=[cosd(theta) -sind(theta) 0; sind(theta) cosd(theta) 0; 0 0 1];
            T=R0*T0*inv(R0);
            imz=stbin(:,:,iz);
            im1=stbin(:,:,i1);
            me1=mean(im1(:));
            im1=im1-me1;
            
            
            mez=mean(imz(:));
            imzc=AffineTransform(imz-mez,T);
            subplot(221);
            imacs(imzc);
            title(num2str([i1 angles(i1)]));
            cc=fftshift(real(ifftn(fftn(im1).*conj(fftn(imzc)))))/nw^2-prior;
            cc=GaussFilt(cc,fcc);
            subplot(222);
            imacs(cc);
            hold on;
            [pk, ix, iy]=max2di(cc);
            plot(ix,iy,'bs');
            plot(pars(:,1)+ctr,pars(:,2)+ctr);
            hold off;
            subplot(223);
            q=circshift(im1,-round([ix-ctr iy-ctr]));
            imacs(q);
            shStack(:,:,i1)=q;
            title(num2str(p));
            drawnow;
            p=Simplex(-pk);
            
        end;
        p=Simplex('centroid');
        par1=[ix-ctr iy-ctr p];
        pars(i1,:)=par1;
        disp([i1 par1]);
    end;
    subplot(224);
    plot(pars(ist:ien,1),pars(ist:ien,2),'.-','markersize',10);
end;

%%  Manual alignment of images
if doManualAlign
    figure(1);
    SetGrayscale;
    subplot(1,1,1);
    ivals2=[iz:-1:1 iz+1:nim];
    pars1=pars;
    
    for i1=ivals2
        imacs(stbin(:,:,i1));
        [ix1,iy1]=ginput(1);
        pars1(i1,1:2)=-[ix1 iy1]+ctr;
    end;
    
    store pars.mat pars pars1 stbin
end;
%% Refine the manual alignment to an ice ball
load pars.mat

for i=1:nim
    shStack(:,:,i)=circshift(stbin(:,:,i),round(pars1(i,1:2)));
end;

imovie(shStack,.2);
%%
% Refine the ice-ball position
ndis=64;
nw=size(stbin,1);
ref=Crop(sum(shStack,3),64);
msk=1-fuzzymask(64,2,18);
imacs(ref+1000*msk);
me=(msk(:)'*ref(:))/(msk(:)'*msk(:));
ref=ref-me;
plot(sect(ref))

% fref=conj(fftn(ref));
fref=conj(fftn(Crop(ref,nw*2)));

sh2Stack=zeros(nw*2,nw*2,nim,'single');
p2=pars1;
for i=1:nim
    %     fq=fftn(Crop(shStack(:,:,i),64));
    fq=fftn(Crop(stbin(:,:,i),2*nw));
    P=FourierShift(nw*2,pars1(i,1:2)-p20(i,1:2));
    cc=fftshift(real(ifftn(P.*fq.*fref)))/(4*nw^2);
    cc=Crop(cc,64);
    [val,p2(i,1),p2(i,2)]=max2di(cc);
    imacs(cc);
    drawnow;
    sh2Stack(:,:,i)=Crop(real(ifftn(P.*fq)),nw*2);
end;
p2(:,1:2)=p2(:,1:2)-33;
%%

shcStack=sh2Stack(525:1244,884-29:1546+28,:);
%
for i=1:nim
    q=shcStack(:,:,i);
    msk=single(q>100);
    me=(msk(:)'*q(:))/sum(msk(:));
    shcStack(:,:,i)=msk.*(q-me);
end;

imovie(shcStack);

%% Try for reconstruction
nst=size(shcStack,1);
nsd=nst/6;
% nsd=nst/4;
dsStack=Downsample(shcStack,nsd,1);
%
fvol=gridMakeNullFT(nsd,3);

degAngles=zeros(nim,3);
degAngles(:,1)=-90;
degAngles(:,2)=angles;
degAngles(:,3)=90;

eAngles=rsDegToEuler(degAngles);
    fvol=gridMakeNullFT(nsd,3);
    
    for i=1:nim
        proj=(dsStack(:,:,i));
        imacs(proj);
        title(num2str([i 180*eAngles(i,1)]));
        drawnow;
        nslice=gridMakePaddedFT(proj);
        fvol=gridInsertPlane(nslice,fvol,eAngles(i,:));
    end
%    %
r=Radius3(size(fvol.PadFT,1));
fs=.2*fvol.np;
filt=r.*exp(-r.^2/(2*fs^2));
fvol1=fvol;
fvol1.PadFT=fvol.PadFT.*filt;

vol=gridRecoverRealImage(fvol1);

ShowSections2(vol);
%%
subplot(339);
imovie(vol(:,:,40:100),.2);

