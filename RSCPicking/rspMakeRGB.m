function cimg=rspMakeRGB(dis,img,masks)
% given the already-scaled img(nx,ny) create the uint8 array cimg(nx,ny,3)
% which has colored mask overlays.  cimg is ready to be shown by image().
% masks has values in the range 0..1 for
% 1: vesicle mask (white vesicles on black background)
% 2: bad vesicle mask
% 3: overall mask (white = forbidden regions)
% The dis structure has fields
% org(x,y)
% size(x,y)
% mode (1..5: orig, subtr, ghost, cc, var)
% ghostColor  (1x3)
% ghostColorBad
%
[nx, ny]=size(img);
if any([nx ny]==0)
    cimg=zeros([dis.size 3],'single');
    return
end;


% Make the gray-scale image
cfim=repmat(rot90(img),[1 1 3]);  % copy into R G B channels

% Add in the vesicle ghosts
if dis.showGhosts==1
    gColor=(1-dis.ghostColor)*dis.ghostAmp; % color to subtract for membrane
    oColor=(1-dis.overlapColor);
    if size(masks,3)<5
        masks(1,1,5)=0;
    end;
    for i=1:3
        cfim(:,:,i)=cfim(:,:,i).*(1-rot90(masks(:,:,1))*gColor(i))...
                                .*(1-rot90(masks(:,:,5)*oColor(i)));
    end;
    bColor=(1-dis.ghostColorBad)*dis.ghostAmp; % color for bad vesicle
    for i=1:3
        cfim(:,:,i)=cfim(:,:,i).*(1-rot90(masks(:,:,2))*bColor(i));
    end;
end;

% Add in masks
cimg=uint8(zeros(ny,nx,3));
if dis.showMask
    mColor=(1-dis.maskColor);
    bColor=(1-dis.blankColor);
    for i=1:3 % loop over colors
        cimg(:,:,i)=uint8(cfim(:,:,i).*(1-rot90(masks(:,:,3)).*mColor(i))...
            .*(1-rot90(masks(:,:,4))*bColor(i)));
    end;
else
    cimg=uint8(cfim);
end;

% crop the display coordinates
dis.org=max(0,dis.org);  % origin is zero-based
dis.org=min([nx ny]-dis.size,dis.org);
% Crop the image to the display size, taking into account the image display
x0=ny-dis.size(2)-dis.org(2)+1; % actually, bottom
x1=ny-dis.org(2); % top
y0=dis.org(1)+1;  % left
y1=dis.org(1)+dis.size(1);  % right
cimg=cimg(x0:x1,y0:y1,:);
