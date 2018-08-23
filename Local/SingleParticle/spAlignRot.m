function [sp aliImgs]=spAlignRot(stack,sp,refs,doFlip)
% function [sp aliImgs]=spAlignRot(stack,sp,refs,doFlip)
% Perform rotational multi-reference alignment of a stack of images, using
% polar coordinates. To the input sp structure are added the fields
% sp.rot    nim x 1 angles, in degrees
% sp.flip   nim x 1 boolean
% sp.class  nim x 1, index of best-match reference
% sp.cc     nim x 1, peak cc value
%  The fast polar-coordinate rotation.

if nargin<4
    doFlip=1;
end;
nflip=1+(doFlip>0);

[n ny nims]=size(stack);
[nr nry nrefs]=size(refs);
display(['spAlignRot: ' num2str(nims) ' images, ' num2str(nrefs) ' references']);
if nr ~= n % need to change the size of the refs
    refs=spDownsample(refs,n);
end;
% normalize the references
for i=1:nrefs
    r=refs(:,:,i);
    refs(:,:,i)=r/sqrt(r(:)'*r(:));
end;
prefs=gridToPolar(refs);
disp('  Converting images to polar coordinates.')
pause(0.01);  % force update
pimgs=gridToPolar(stack);
ntheta=size(prefs,2);


disp('  Rotational correlation');
pause(.01);

mx2=zeros(nims,1);
    j2=mx2;
    ang2=mx2;
    k2=mx2;

parfor i=1:nims
    mx0=zeros(nflip,1);
    ang0=mx0;
    mx1=zeros(nrefs,1);
    ang1=mx1;
    k1=mx1;
    ccs=spRotCorr(pimgs(:,:,i),prefs);
    for j=1:nrefs
        for k=1:nflip
            cc=ccs(:,1,j,k);
            cc(ntheta+1:ntheta+2)=cc(1:2);
            [mx0(k) ang0(k)]=max1di(cc);  %maximize over angles
        end;
        [mx1(j) kv]=max(mx0);  % max over flipping
        ang=ang0(kv);
        if kv>1
            ang=1.5*ntheta-ang+2;
        end;
        ang1(j)=ang;
        k1(j)=kv;
    end;
    [mx2(i) j2(i)]=max(mx1);  % maximize over reference
    ang2(i)=ang1(j2(i));
    ang2(i)=mod(ang2(i),ntheta);
    k2(i)=k1(j2(i));
    if (mod(i,1000)==0)
        disp(i);
        pause(.01);
    end;
end;

sp.rot=(ang2-1)*360/ntheta;
sp.flip=k2-1;
sp.class=j2;
sp.cc=mx2;

if nargout>1 % if the user requested the aligned stack
    disp('  Final rotation of images');
    pause(.01);
    %     aliImgs=TransformImages(stack,xytfr);
    [sp0 aliImgs]=spTransformStack(stack,sp);
end;
disp('done.');
