function m=meMakeModelVesicleDiscs(mi,n,vindex,risetime)

% function v=meMakeModelVesicles(mi,n,vindex,risetime)
% function v=meMakeModelVesicles(mi,scl,vindex,risetime)
% Given the info structure mi, make a model consisting of discs at each
% vesicle position
% for each vindex value in the mi.vesicle arrays; default is every one for
% which mi.vesicle.ok(index,1)==true.
% The result is an n-sized image (may be rectangular)
% computed as a possibly downsampled version of the original micrograph.
% If the second argument is a struct we assum it's
%   scl.n size of the output image
%   scl.M image scale matrix from meGetImageScale() which is used in
%   interpreting the vesicle coordinates.

if isa(n,'struct')% has n and an affine matrix
%     dsShift=-n.M(1:2,3)';
    M=n.M;
    ds=M(1,1);
    n=n.n;
else % in th other cases, we assume mi.imageSize is simply a multiple of n
%     And we construct M with no shift.
    ds=mi.padImageSize(1)/n(1);   % downsample factor
    M=[ds 0 0; 0 ds 0; 0 0 1]; 
end;

v=zeros(n,'single');  % default is a zero image.
if numel(mi.vesicle.s)<1
    return
end;
badS=isnan(mi.vesicle.s(1));
mi.vesicle.ok(badS,:)=false;

% Get the membrane cross-section density.
% If no model is present, return a zero image.

nv=numel(mi.vesicle.x);
if nv<1
    return
end;

if nargin<3
    vindex=0;
end;
if numel(vindex)<1
    return
end;
if vindex(1)==0
    vindex=find(mi.vesicle.ok(:,1));  % all found vesicles.
end;
vindex = vindex(vindex<=nv);  % don't allow out-of-range indices

if numel(n)~=2
    n=[1 1]*n;
end;

nim=numel(vindex);
sumv=single(zeros(n));

% Transform from original micrograph coordinates to local coords
globalXY=[mi.vesicle.x mi.vesicle.y ones(nv,1)]';

vesXY=M\globalXY; % zero-based local coordinates.

for k=1:nim
    ind=vindex(k);
    % Get the coordinates and radius, scaled down by ds    
        vx=vesXY(1,ind);
        vy=vesXY(2,ind);
vr=mi.vesicle.r(ind,:)/ds;
t=risetime/ds;

    % Accumulate the density
    sumv=sumv+meMakeVesicleDisc(n,vr,t,[vx vy]);
end;
m=sumv;
