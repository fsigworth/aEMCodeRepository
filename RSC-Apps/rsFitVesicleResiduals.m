% rsFitVesicleResiduals
% Experiment with fitting vesicle residuals in a particle stack.  A
% classification without alignment is made of the original stack; then the
% residual-subtracted stack is averaged using the same classification to
% allow direct comparison.
writeOutputFile=1;  % Write the *stackc.mrc output file
usePhysicsModel=0;
particleWidth=100;  % angstroms
maxHalfAngle=60;  % degrees
maxHalfWidth=200; % angstroms

membraneWidth=60;

[fname, pa]=uigetfile('*si.mat','Select a stack info file *si.mat');
if isnumeric(pa) % File selection cancelled
    return
end;
cd(pa);

load(fname);
p=strfind(fname,'si.mat');
baseName=fname(1:p-1);
imgs=ReadMRC([baseName 'stack.mrc']);
imgsc=zeros(size(imgs),'single');  % corrected image stack
fitsc=zeros(size(imgs),'single');
% imgs=ReadMRC([baseName 'ustack.mrc']);
nim=numel(si.miIndex);
if nim ~= size(imgs,3)
    error(['Inconsistent number of images ', num2str([nim size(imgs,3)])]);
end;
numImages=nim
n0=size(imgs,1);  % raw image size
n1=NextNiceNumber(n0*1.1);  % pad to 10% bigger
% [r th]=Radius(n1);
% theta=90-th*180/pi;

ctr=n1/2+1;

ds=si.pixA / si.mi{si.miIndex(1)}.pixA;  % downsampling factor
% multiply the vesicle model by the voxel size.

pw=particleWidth/si.pixA;
mw=membraneWidth/si.pixA;

% v=single(real(ifftn(fftn(double(sumv)).*ifftshift(double(H)))));  % Filter with the ctf.
%%
v=zeros(n1,n1,nim);
sum=zeros(n1,n1);
sum0=zeros(n1,n1);

if usePhysicsModel
    nrefs=4;  % derivative basis functions
else
    nrefs=12;
end;
coeffs=zeros(nrefs,nim);  % shell basis functions
oldimi=0;
mask0=Crop(SquareWindow(n0),n1);
figure(1);
SetGrayscale;
for i=1:nim
    imi=si.miIndex(i);
    if imi ~= oldimi
        mi0=si.mi{imi};
        H=meGetEffectiveCTF(mi0,n1,ds);
        iH=ifftshift(H.*meGetNoiseWhiteningFilter(mi0,n1,ds)); % zero-origin filter
        vd=meDownsampleVesicleModel(mi0.vesicleModel,ds)*ds*mi0.pixA;
        oldimi=imi;
    end;
    coords=mi0.particle.picks(si.miParticle(i),:);
    iVesicle=coords(4);
    r=si.rVesicle(i);
    s=si.sVesicle(i);
    y0=si.yClick(i);

    [rs, th]=Radius(n1,[ctr ctr-y0]);
    theta=90-th*180/pi;
    w=5;
    innerAngle=25;
    outerAngle=80;
    mask=0.5*(erf((abs(theta)-innerAngle)/w)-erf((abs(theta)-outerAngle)/w));
    mask=mask.*mask0;
    img=Crop(imgs(:,:,i),n1);
    subplot(2,3,1); imacs(img);
    title([num2str(i) ' / ' num2str(nim)]);
%     subplot(2,3,3); imacs(img.*mask);
%         subplot(2,3,4); imacs(mask); 
        drawnow;
    mski=mask(:);
    
if usePhysicsModel    
%%      Make the 3 vesicle-derivative basis functions
        t2=.05/sqrt(r);
    dr=.25;
        v1=s*VesicleFromModel(n1,r,vd,[ctr ctr-y0],0);
        v2=s/dr*VesicleFromModel(n1,r-dr,vd,[ctr ctr-y0],0);
        v4=s/dr*VesicleFromModel(n1,r+dr,vd,[ctr ctr-y0],0);
        v3=s*VesicleFromModel(n1,r,vd,[ctr ctr-y0],t2);
        v0=0*v1+1;
        v1f=real(ifftn(fftn(v1).*iH));
        v2f=real(ifftn(fftn(v4-v2).*iH)); % first derivative
        v3f=real(ifftn(fftn(v3-v1).*iH))/t2; % curvature change
        f=[v0(:) v1f(:) v2f(:) v3f(:)];
        fMasked=double([mski.*v0(:) mski.*v1f(:) mski.*v2f(:) mski.*v3f(:)]);
    subplot(2,3,5); imacs(v3f.*mask);
% % img=v3f;
else
% Alternatively make the 12 shell basis functions
    v=zeros(n1,n1,nrefs);
    dx=mw/nrefs;
    nx=ceil(dx*5)*2+1;  % 5 stds; odd number of points
    temp0=zeros(nx,1);
    temp0((nx+1)/2)=1/nx;
    template=GaussFilt(temp0,.2/dx);
    for j=1:nrefs
        dr=(j-(nrefs-1)/2)*dx;
        v0=s*VesicleFromModel(n1,r+dr,template,[ctr ctr-y0],0);
        v(:,:,j)=real(ifftn(fftn(v0).*iH));
    end;
%%    
    f=reshape(v,n1*n1,nrefs);   
    fMasked=double(f.*repmat(mski,1,nrefs));
end;
    y=double(mski.*img(:));
    
    a=LinLeastSquares(fMasked,y);
    coeffs(:,i)=a;
    a1=a;
    fit=reshape(f*a1,n1,n1);
    
    sum=sum+img-fit;
    subplot(2,3,4); imacs(sum);
    title('Sum(corrected)');
    
    sum0=sum0+img;
    subplot(2,3,3); imacs(sum0);
    title('Sum(original)');
    
    subplot(2,3,2); plot(a);
%     axis([0 inf -50 50]);
    drawnow;
    imgsc(:,:,i)=Crop(img-fit,n0);
    fitsc(:,:,i)=Crop(fit,n0);
end;

%%
%
% Initial clssification
mskRadius=0.3*n0;
rmask=fuzzymask(n0,2,mskRadius,mskRadius/5);
figure(1);
subplot(2,3,1);
imacs(rmask);
title('MSA mask');

nfactors=20;
nclasses=100;
%  Classifying the original images
[means0, nm, inds, fvecs]=Classifier(imgs,nclasses,nfactors,2,rmask);  % 20 factors, decimate by 2.
%
figure(4);
ImagicDisplay1(means0);  % original classes
%%

% Compute full-sized means
means1=single(zeros(n0,n0,nclasses));
imgsc1=imgsc-0*fitsc;  % corrected images, with over-correction
for i=1:nclasses
    means1(:,:,i)=squeeze(mean(imgsc1(:,:,inds==i),3));
end;
figure(5);
ImagicDisplay1(means1);  % corrected classes

meansf=single(zeros(n0,n0,nclasses));
for i=1:nclasses
    meansf(:,:,i)=squeeze(mean(fitsc(:,:,inds==i),3));
end;
figure(6);
ImagicDisplay1(meansf);  % residuals subtracted


%%
if writeOutputFile
    WriteMRC(imgsc,si.pixA,[baseName 'stackc.mrc']);
end;
