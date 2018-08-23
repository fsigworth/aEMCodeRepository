function ri=re2DMakeRunInfoScaled(ri,n)
% function ri=reMakeRunInfoScaled(ri,n)
% Make pixel-size-dependent entries into the ri structure, based on
% the current image size n.

n=single(n);
ri.nCurrent=n;
ds=ri.nCropU/n;  % downsampling factor relative to the original stack
ri.pixA=ri.pixAU*ds;

% % Derived angle parameters: downsample only the alphas
riTemp=reSetRefAngles(ri.angleStepsU*ds,ri.angleLimits,ri.isos,false,ri);
ri.alphas=riTemp.alphas;
%ri.alphasI=riTemp.alphasI;
ri.angleN(1)=riTemp.angleN(1);
ri.angleStep(1)=riTemp.angleStep(1);

% Masks
volMaskRadius=.25*n;
volMaskHeight=round(.35*n);
volMask=zeros(n,n,n,'single');
ctr=floor(n/2+1);
for iz=ctr-volMaskHeight:ctr+volMaskHeight
    volMask(:,:,iz)=fuzzymask(n,2,volMaskRadius,.5);
end;

ri.volMask=GaussFilt(volMask,.02*ds);
ri.radius=floor(n/2)-ri.maxShiftU/ds;
nt=2*ceil(ri.maxShiftU/ds)+1;  % number of translation steps in x or y
ri.radius=floor((n-nt)/2);
ri.nTrans=nt;
ri.softMask=fuzzymask(n,2,ri.radius*.9,ri.radius*.2);

fscMaskRadius=.33*n;
fscMaskWidth=.05*n;
ri.fscMask=fuzzymask(n,3,fscMaskRadius,fscMaskWidth);
