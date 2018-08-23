function v=meMakeModelVesicles(mi,n,vindex,doCTF,doPW)
% function v=meMakeModelVesicles(mi,n,vindex,doCTF,doPW)
% function v=meMakeModelVesicles(mi,img,vindex,doCTF,doPW)
% Given the info structure mi, make a CTF-filtered, scaled vesicle model
% for each vindex value in the mi.vesicle arrays; default is every one for
% which mi.vesicle.ok==1.
% The result is an n-sized image (may be rectangular)
% computed as a possibly downsampled version of the original micrograph.
% By default the CTF is applied (doCTF=1) but the prewhitening from the
% noise model is not (doPW=0).
% If vindex is not given or is zero, vesicles are created only for rows of
% mi.vesicle.ok that are all 1s.
% If the second argument is an image, then n is taken as size(img) and the
% final v is scaled to match img by least-squares.


if numel(n)>2 % we supplied an image
    m0=n;
    n=size(n);
    doFitImage=1;  % We'll do least-squares to the image.
else
    doFitImage=0;
end;

v=single(zeros(n));  % default is a zero image.
ds=mi.imageSize(1)/n(1);   % downsample factor
if numel(mi.vesicle.s)<1
    return
end;
badS=isnan(mi.vesicle.s(1));
mi.vesicle.s(badS,:)=0;

v=single(zeros(n));  % default, return zeros.

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
if vindex(1)==0 || ~isfield(mi.vesicle,'ok')
    mi.vesicle.ok=true(numel(mi.vesicle.x),4);
    mi.vesicle.ok(badS,:)=false;
    vindex=find(all(mi.vesicle.ok,2));  % ok is an n x 4 or so matrix.
end;
vindex = vindex(vindex<=nv);  % don't allow out-of-range indices

if nargin<4
    doCTF=1;
end;
if nargin<5
    doPW=0;
end;

if numel(n)~=2
    n=[1 1]*n;
end;

if ~(isfield(mi,'vesicleModel') && numel(mi.vesicleModel)>1)
    error('No vesicle model');
end;

% multiply the vesicle model by the voxel size.
vd=meDownsampleVesicleModel(mi.vesicleModel,ds)*ds*mi.pixA;

nim=numel(vindex);
sumv=single(zeros(n));
for k=1:nim
    ind=vindex(k);
    % Get the coordinates and radius, scaled down by ds
    vx=(mi.vesicle.x(ind))/ds+1;  % zero-based coordinate
    vy=(mi.vesicle.y(ind))/ds+1;
    vr=mi.vesicle.r(ind,:)/ds;
    
    % Accumulate the vesicle density
%     sumv=sumv-mi.vesicle.s(ind)*VesicleFromModel(n,vr,vd,[vx vy]);
    sumv=sumv-VesicleFromModelGeneral(n,vr,vd,[vx vy],mi.vesicle.s(ind,:,1));

    % -------------Create extra ring components--------------
if isfield(mi.vesicle,'extraPeaks') && numel(mi.vesicle.extraPeaks)>0
    np=numel(mi.vesicle.extraPeaks);
    if size(mi.vesicle.s,3)==np+1 % extra amps are part of the s array
        exAmps=ds*mi.pixA*shiftdim(mi.vesicle.s(ind,:,2:end),1);
    elseif isfield(mi.vesicle,'extraS') && size(mi.vesicle.extraS,3)==np+1  
            % old version, use the extraS field
        exAmps=ds*mi.pixA*shiftdim(mi.vesicle.extraS(ind,:,:),1);
    else
        exAmps=0;
    end;    
    if any(exAmps(:))
        exPos=mi.vesicle.extraPeaks/ds;
        exSigma=mi.vesicle.extraSD/ds;
        sumv=sumv-VesicleFromRings(n,exPos,exSigma,vr,[vx vy],exAmps);
    end;
end;

end;
v=sumv;  % default returned value
H=1;
if doCTF % operate with the CTF
    H=meGetEffectiveCTF(mi,n,ds);
end;
if doPW  % Pre-whitening filter
    H=H.*meGetNoiseWhiteningFilter(mi,n);
end;
if doCTF || doPW  % do the filtering
    v=single(real(ifftn(fftn(double(sumv)).*ifftshift(double(H)))));  % Filter with the ctf.
end;

if doFitImage % do least-squares fit to supplied image
    F=[ones(numel(m0),1) v(:)];
    a=LinLeastSquares(F,m0(:));
    v=v*a(2);
end;



