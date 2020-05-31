function [sclOut,rgb]=imacx2(m,power,pars)
% function scl=imacx2(m,power,pars)
% Display a complex-valued image, using Lab color.
% Special case: in pixels where the imaginary part is inf, a grayscale
% value is derived from the real part, in the range 0..1.
% pars is a struct containing extra parameters. 
%  pars.scl --scaling of m, normally chosen so max(abs(m*scl))==1
%  pars.sat --color saturation. 1 is normal.
%  pars.antiAliasing --If nonzero, the brightness and colors are each downsampled
%    to reduce aliasing in the display. Values >1 represent fractions of
%    Nyquist cutoff in the Fourier downsampling.
%  pars.x
%  pars.y
% function scl=imacx2(x,y,m,power,scl,sat,doAntiAliasing)
% function scl=imacx2(m,power,scl,sat,doAntiAliasing)
% function scl=imacx2(x,y,m)
if nargin<3
    pars=struct;
end;
if nargin<2
    power=1;
end;


dpars.scl=[];
dpars.sat=1;
dpars.antiAliasing=0;
dpars.x=[];
dpars.y=[];

pars=SetDefaultValues(dpars,pars);

% Autoscaled plot of complex-valued matrix m.  Intensity of the plot is
% determined by abs(m), and the color represents the phase.
% If the argument 'power' is given, then the intensity is abs(m).^power.
% Typically power<1 to expand the dynamic range of the display.
% If desired, the magnitude scaling is returned, for use with further
% displays.  The scaling is computed as scl=1/max(abs(m(:)).^power).
% x and y arguments, if given, are 1D vectors specifying the extent of the
% plotted points.
%
% [x,y]=ndgrid(1:1023);
% lambda=16; % wavelength in pixels
% m=(exp(1i*2*pi*y/lambda));
% m=m*1;
% 

m=single(squeeze(m))';
sz=size(m);
monoMask=isinf(m);
r1=ones(sz,'single');
r1(~monoMask)=abs(m(~monoMask));
scl=pars.scl;
if numel(scl)<1
    scl=1/max(r1(:));
end;

if nargin<2 || numel(power)<1
    pars.power=1;
end;

re1=real(m)./r1;
% re1(monoMask)=real(m(monoMask));
ri1=imag(m)./r1;
ri1(monoMask)=0;
ar1=(scl*r1).^power;
if pars.antiAliasing
    h=gca;
    h.Units='pixels';
    plotSize=round(h.Position(4:-1:3));
    h.Units='normalized';
%     fc=plotSize./(pars.antiAliasing*sz);  % 2D cutoff frequency
    fw=plotSize/(3*pars.antiAliasing);
    H=ifftshift(fuzzymask(sz,2,fw,fw/5));
    re1=real(ifftn(fftn(re1).*H));
%     ri1(monoMask)=0;
    ri1=real(ifftn(fftn(ri1).*H));
    ar1=real(ifftn(fftn(ar1).*H));

%     ar1=Downsample(ar1,plotSize,0,msk);
%     re1=Downsample(re1,plotSize,0,msk);
%     ri1=Downsample(ri1,plotSize,0,msk);
%     ar1=Downsample(ar1,plotSize,0,msk);
% else
%     plotSize=sz;
end;

rLimit=70;
aLimit=40*pars.sat;
bLimit=40*pars.sat;
rgb=lab2rgb([rLimit*ar1(:) aLimit*re1(:) bLimit*ri1(:)]);
rgb(monoMask(:),:)=repmat(real(m(monoMask(:))),1,3);
rgb=reshape(rgb,[sz 3]);

% %     'ColorSpace','linear-rgb');
% for k=1:3
%     rgb(monoMask,k)=real(m(monoMask));
% end;
if numel(pars.x)>0 && numel(pars.y)>0
    image(pars.x,pars.y,rgb);
else
    image(rgb);
end;
axis xy;

if nargout>0
    sclOut=scl;
end;
% 
% % wrap t about zero, not -pi
% % td=t/(2*pi)+(t<0);
% td=t/(2*pi)+.5;
% d=64*(1-eps); % magic factor
% md=257+floor(r*d)+64*floor(d*td);
% 
% % draw the image
% image(md');
% axis xy
