function mColor=ShowImageAndBoxes(m,ctrs,boxSize,boxLine,color,thresh)
% Display a grayscale image overlaid with boxes having the interior
% dimension boxSize and linewidth boxLine.  Color is a vector e.g.
% yellow = [1 1 0].
% ctrs is an nboxes x 2 array giving x,y pixel positions.
% This can also be called as ShowImageAndBoxes(m) in which case it will
% show a grayscale image.
newVersion=1;
newBoxLine=1;

if nargin<6
    thresh=1e-4;
end;
areBoxes=nargin>1 && numel(ctrs)>1;
n=size(m);


if newVersion && nargout<1  % simpler code, using hold command
    mScaled=single(imscale(m,256,thresh));  % scale image
    imaga(mScaled);
    if areBoxes
        hold on;
        nb=size(ctrs,1);
        xvec=[-1 -1 1 1  -1]*boxSize/2;
        yvec=[-1  1 1 -1 -1]*boxSize/2;
        for i=1:nb
            plot(xvec+ctrs(i,1),yvec+ctrs(i,2),'linewidth',newBoxLine,'color',color);
        end;
        hold off;
    end;
    
else  % old version, set pixel values.

    mScaled=single(imscale(m,1,thresh));  % scale to 0..1

    if areBoxes
        if nargin<5
            color=[1 1 0];  % yellow
        end;
        if nargin<4
            boxLine=round(n(1)/512);
        end;
        innerBox=zeros(boxSize);
        xBoxSize=boxSize+2*boxLine;
        box=Crop(innerBox,xBoxSize,0,1);  % fill around it with ones
        blank=zeros(xBoxSize);
        blankMask=1-Mask(zeros(n),ctrs',blank,box);
        drawMask=Mask(zeros(n),ctrs',blank,box);  % opaque boxes
        % drawMask=Mask(zeros(n),ctrs,white,box);  % overlapping boxes
    end;
    
    mColor=single(zeros([n 3]));
    if areBoxes
        for i=1:3
            mColor(:,:,i)=mScaled.*blankMask+color(i)*drawMask;
        end;
    else
        mColor=repmat(mScaled,1,1,3);
    end;
    
    if nargout<1
        mColorR=zeros(n(2),n(1),'single');
        for i=1:3
            mColorR(:,:,i)=rot90(mColor(:,:,i));
        end;
        image(mColorR);
    end;
    
end;