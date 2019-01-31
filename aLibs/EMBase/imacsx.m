function scale=imacsx(m,power,scl)
%  imacx2(m), imacsx(m,power,scl)
% Autoscaled plot of complex-valued matrix m.  Intensity of the plot is
% determined by abs(m), and the color represents the phase.
% If the argument 'power' is given, then the intensity is abs(m).^power.
% Typically power<1 to expand the dynamic range of the display.
% If desired, the magnitude scaling is returned, for use with further
% displays.  The scaling is computed as scl=1/max(abs(m(:)).^power).

m=single(squeeze(m')); % take the transpose for Cartesian coords
%  make the argument single for image command below.

if nargin<2
    r=abs(m);
else
    r=abs(m).^power;
end;

t=mod(atan2(imag(m),real(m))/(2*pi),1);  % angle from 0 to 1-eps
mx=max(max(r));
if mx > 0 && nargin<3
    r=r/mx;
    scl=1/mx;
else
    r=min(r*scl,1); % clip if the scale is too large.
end;

s=1-.5*r; % increase brightness at extremes.

image(hsv2rgb(t,s,r));
% 
% % wrap t about zero, not -pi
% % td=t/(2*pi)+(t<0);
% td=t/(2*pi)+.5;
% d=64*(1-eps); % magic factor
% md=257+floor(r*d)+64*floor(d*td);
% 
% % draw the image
% image(md');
axis xy
if nargout>0
    scale=scl;
end;
