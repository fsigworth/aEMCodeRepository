% MakeMosaic.m

ds=8;  % binning factor
refIntensity=2100;
iceConstant=1;  % e-fold decay for 1um ice
thicknessRange=.4;  % um of ice

% ----------------Get the stack files---------------

[fname, pa]=uigetfile('*.st.mdoc','Select a mosaic doc file');
if isnumeric(pa)  % user clicked Cancel
    return
end;

cd(pa);

%%
stackName=fname(1:end-5);
disp(['Reading ' stackName]);
[st, pixA]=ReadEMFile(stackName);
[nx, ny, nz]=size(st);
nxw=floor(nx/ds);
nyw=floor(ny/ds);
disp('Downsampling');
% the image size had better be a multiple of ds
fst=Downsample(Crop(single(st),[nxw nyw]*ds,1),[nxw nyw],1);

mid=mean(fst(:));


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
coords=zeros(0,2);
while ~feof(mdoc)
    if numel(s)>0
        c=textscan(s,'%s = %f %f %f','emptyvalue',0);
        tok=char(c{1});
        switch tok
            case '[ZValue'
                iz=c{2}+1;
            case 'PieceCoordinates'
                xpos=floor(c{2}/ds);
                ypos=floor(c{3}/ds);
                coords(iz,:)=[c{2} c{3}];
                mosaic(xpos+1:xpos+nxw,ypos+1:ypos+nyw)=fst(:,:,iz);
                assigned(xpos+1:xpos+nxw,ypos+1:ypos+nyw)=true;
        end;
    end;
    s=fgetl(mdoc);
end;
mosaic(~assigned)=refIntensity*exp(-thicknessRange*0.7/iceConstant);

% norm. image intensity = exp(-thk/iceConstant)
% expand to 256 when norm intensity=0, 0 when
% normIntensity=exp(-thickRange/iceConstant).
% hence relthk=(-iceConstant*log(intensity))/thicknessRange;
mosaic=GaussFiltDCT(mosaic,.05);


b='d';
displayMode=0;
while (b~='q') && (b~='Q') % q = quit; Q = quit but return to this image later
    switch b
        
        
        case 'd'
            switch displayMode
                case 0
                    scaledMosaic=1+iceConstant*log(max(1,mosaic)/refIntensity)/thicknessRange;
                    xs=(1:nxw)*ds*pixelSpacing/1000;  % coord in nm
                    ys=(1:nyw)*ds*pixelSpacing/1000;  % coord in nm
                    imac(xs,ys,scaledMosaic*200)
                    colormap jet
                    colorbar('yticklabel',1000*(4:-1:0)*thicknessRange/5);
                    xlabel('position, um');
                    b=2;
                case 1
                    SetGrayscale;
                    imacs(xs,ys,mosaic);
                    colorbar;
            end;
            displayMode=mod(displayMode+1,2);
    end;
    
    [cx,cy,b]=ginput(1);
end;
disp('Done.');

return








figure(1);
SetGrayscale;
subplot(2,2,1);
plot(coords(:,1),coords(:,2),'.-','markersize', 10);

[nx, ny, nz]=size(fst);

npad=[2*nx 2*ny];
i=2; j=1;
fstx=Crop(fst,npad,1);  % pad the whole stack

cc=fftshift(real(ifftn(fftn(Crop(fst(:,:,i),npad)).*conj(fftn(Crop(fst(:,:,j),npad))))));

subplot(2,2,2);
imacs(cc);

subplot(2,2,3);
imacs(fstx(:,:,1));
subplot(2,2,4);
imacs(fstx(:,:,2));
plot(sect(fstx(:,:,2)))

return
%%
thresh=700;
padImgs=zeros([npad nz],'single');
padMsks=zeros([npad nz],'single');
for i=1:nz
    m=fst(:,:,i);
    msk=m>thresh;
    %      msk=GaussFilt(msk,.05)>.99;
    mmsk=m.*msk;
    me=sum(mmsk(:))/sum(msk(:));
    mcorr=(m-me).*msk;
    subplot(2,2,3); imacs(m);
    subplot(2,2,4); imacs(mcorr);
    drawnow;
    padImgs(:,:,i)=Crop(mcorr,npad);
    padMsks(:,:,i)=Crop(msk,npad);
end;
fPadImgs=fft2(padImgs);
fPadMsks=fft2(padMsks);
%%
i=1;
j=2;
rotMat=[0 -1; 1 0];
vector=(coords(j,:)-coords(i,:))*rotMat;
[shifts,overlap]=MosaicCC(padImgs,fPadImgs,fPadMsks,i,j);
angle=180/pi*(atan2(vector(2),vector(1))-atan2(shifts(2),shifts(1)));
umPixels=hypot(shifts(2),shifts(1))/hypot(vector(2),vector(1));
predMat=umPixels*[cosd(angle) -sind(angle); sind(angle) cosd(angle)];

vector*predMat
pixCoords=(coords*rotMat)*predMat;
pixCoords=pixCoords-repmat(pixCoords(1,:),nz,1)
return;

nomShift=[930 0; 0 930];
stageShift=[-6.4 -.9; .9 -6.4];
diffCoords=-diff(coords);
relCoords=[coords(:,1)-coords(1,1) coords(:,2)-coords(1,2)];
intRelCoords=round(relCoords*stageShift/(stageShift(:,1)'*stageShift(:,1)));

neighbors=false(nz,nz);
ls=true(nz,nz);
for i=1:nz
    ref=repmat(intRelCoords(i,:),nz,1);
    difs=abs(intRelCoords-ref);
    difs(1:i,:)=0;  % force lower diagonal
    %     neighbors(:,i)=all(difs<2,2) & any(difs,2);
    neighbors(:,i)=sum(difs,2)<2 & any(difs,2);
    
end;
inds=find(neighbors);
[is,js]=ind2sub(size(neighbors),inds);  % js change slowly
ctr=[nx ny]+1;
xds=zeros(nz,nz);
yds=zeros(nz,nz);
for k=1:numel(is);
    i=is(k);
    j=js(k);
    %         compShift=((intRelCoords(j,:)-intRelCoords(i,:))*nomShift)*rotMat;
    [shifts,overlap]=MosaicCC(fPadImgs,fPadMsks,i,j);
    
end;


%         cc0=fftshift(real(ifftn(fPadImgs(:,:,i).*conj(fPadImgs(:,:,j)))));
%         ccm=fftshift(real(ifftn(fPadMsks(:,:,i).*conj(fPadMsks(:,:,j
%         % ccq=cc;
% %         cc=GaussHP(ccq,.02).*fuzzymask(npad,2,100,100,compShift+ctr);
%         cc=GaussHP(ccq,.1);
%         subplot(2,2,2);
%         imacs(cc);
%         [val, ix, iy]=max2d(cc);
%         % cc1=cc;
%         % cc1(ix-10:ix+10,iy-10:iy+10)=0;
%         % val1=max2d(cc1);
%         % ij=[i j]
%         % peaks=[val val1]
%         overlap=ccm(ix,iy);
%         shift=[ix iy]-ctr;
%         disp([j i round(shift/10) round(overlap/1e4);
%             j i round(compShift/10) round(val)]);
%         disp(' ');
%         %  shift=[1180 1500]-ctr;
%         % shift=[-8 912]
% %         shift
%         im2=circshift(padImgs(:,:,j),shift);
%         subplot(2,2,3);
%         im1=padImgs(:,:,i);
%         imacs(im1);
%         subplot(2,2,4);
%         imacs(im2);
%         subplot(2,2,1);
%         imacs(im2+im1);
%         title([i j]);
%         drawnow;
