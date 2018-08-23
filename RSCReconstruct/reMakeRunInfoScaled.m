function ri=reMakeRunInfoScaled(ri,n)
% function ri=reMakeRunInfoScaled(ri,n)
% Make pixel-size-dependent entries into the ri structure, based on
% the current image size n.
if nargin<2
    n=ri.nCurrent;
else
    n=single(n);
    ri.nCurrent=n;
end;
ds=ri.nCropU/n;  % downsampling factor relative to the original stack
ri.pixA=ri.pixAU*ds;

% Scale the timeout values as n^2.  baseTimeout is typically 30 min.
ri.timeout=ri.baseTimeout*(.2+(n/64)^2+(n/64)^3);

% Derived angle parameters
if ri.flags.mode2D
    angleSteps=ri.angleStepsU;
    angleSteps(1)=ri.angleStepsU(1)*ds;
else
    angleSteps=ri.angleStepsU*ds;
end;
ri=reSetRefAngles(angleSteps,ri.angleLimits,ri.isos,false,ri);
ri.angles=reGetAngleList(ri,false);

% Masks
ctr=floor(n/2+1);

% 3D mask
switch ri.volMaskType
    case 'cylinder'
        % Cylindrical mask
        volMaskRadius=.25*n;
        volMaskHeight=round(.35*n);
        volMask=zeros(n,n,n,'single');
        for iz=ctr-volMaskHeight:ctr+volMaskHeight
            volMask(:,:,iz)=fuzzymask(n,2,volMaskRadius,.5);
        end;
        ri.volMask=GaussFilt(volMask,.02*ds);
        
    case 'sphere+ellipsoid'  %Sphere+ellipsoid mask
        ri.volMask=min(1,fuzzymask(n,3,[.32 .32 .15]*n,.1*n,[ctr ctr ctr+.18*n])...
            +fuzzymask(n,3,0.22*n,0.1*n,[ctr ctr ctr-.06*n]));
    case 'sphere'
        ri.volMask=fuzzymask(n,3,0.22*n,0.1*n,[ctr ctr ctr-.06*n]);
    otherwise  % sphere+ellipsoid anyway
        disp(['Volume mask type not recognized: ' ri.volMaskType]);
        ri.volMask=min(1,fuzzymask(n,3,[.32 .32 .15]*n,.1*n,[ctr ctr ctr+.18*n])...
            +fuzzymask(n,3,0.22*n,0.1*n,[ctr ctr ctr-.06*n]));
end;

% 2D masks
ri.radius=floor(n/2)-ri.maxShiftU/ds;
nt=2*ceil(ri.maxShiftU/ds)+1;  % number of translation steps in x or y
ri.radius=floor((n-nt)/2);
ri.nTrans=nt;
ri.softMask=fuzzymask(n,2,ri.radius*.9,ri.radius*.2);

fscMaskRadius=.33*n;
fscMaskWidth=.05*n;
ri.fscMask=fuzzymask(n,3,fscMaskRadius,fscMaskWidth);
