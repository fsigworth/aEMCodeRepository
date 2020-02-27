% rsFitVesicleDeformation(mi,msub)
% vindex=4;
% diffIm0=diffIm;
% diffIm=diffIm0;
% diffIm=real(z.^3.*vfilt1)+randn(size(vfilt20))*.002;

vindex=ind;

ns=size(ms,1);  % downsampled micrograph
ndis=72;

ds=mi.imageSize(1)/ns;

vx=(mi.vesicle.x(vindex))/ds+1;  % assume zero-based coordinates
vy=(mi.vesicle.y(vindex))/ds+1;
vr=mi.vesicle.r(vindex)/ds;
vs=mi.vesicle.s(vindex);

% Get the membrane cross-section density
vd1=mi.vesicleModel;
vd=meDownsampleVesicleModel(vd1,ds)*mi.pixA*ds;
approxShift=round([vx vy]);  % center of cropped image in ms
dShift=[vx vy]-approxShift+ndis/2;
% this is the approximate position of the vesicle relative to
% the center of the full-size image

% Get the filter function
H=ifftshift(meGetEffectiveCTF(mi,ndis,ds));

% original vesicle
v0=-vs*VesicleFromModel(ndis,vr,vd,dShift);
vfilt0=real(ifftn(fftn(v0).*H));
vfiltNorm=vfilt0(:)'*vfilt0(:);
vfilt0=vfilt0/sqrt(vfiltNorm);  % NCC
subplot(223);
imacs(vfilt0);

% differential vesicle
% difference in image if the vesicle radius were 1 original pixel larger.
dr=1/4;  % step of differential 
v1=-vs/(dr*ds)*(VesicleFromModel(ndis,vr+dr,vd,dShift)...
    -VesicleFromModel(ndis,vr-dr,vd,dShift));
vfilt1=real(ifftn(fftn(v1).*H));
vf1Cross=vfilt1(:)'*vfilt0(:);
vfilt1=vfilt1-vf1Cross*vfilt0;  % orthogonalize
vf1Norm=vfilt1(:)'*vfilt1(:);
vfilt1=vfilt1/sqrt(vf1Norm);
subplot(224);
imacs(vfilt1);

% Get inner products
% d0=diffIm(:)'*vfilt0(:);
% d1=diffIm(:)'*vfilt1(:);

% Get functions of theta

[r theta]=Radius(ndis);
z=exp(1i*theta);
plot([real(sect(z)) imag(sect(z))]);
%%
z1=z(:);
nterms=6;
q=complex(zeros(nterms,2));
for i=1:nterms
    q(i,1)=(1+(i>1))*vfilt0(:)'*(z1.^-(i-1).*diffIm(:));
    q(i,2)=(1+(i>1))*vfilt1(:)'*(z1.^-(i-1).*diffIm(:));
end;

% Try making a better model
subplot(221);
imacs(diffIm);

vfilt2=0*vfilt0;
for i=1:nterms
    vfilt2=vfilt2+real(z.^(i-1).*(vfilt0*q(i,1)+vfilt1*q(i,2)));
end;
subplot(224);
imacs(vfilt2);
subplot(222);
imacs(diffIm-vfilt2);
% plot([sect(diffIm) sect(diffIm-vfilt2)])
vindex
vr
abs(q)
subplot(221)
title(round(sum(abs(q(:)))/vs));

