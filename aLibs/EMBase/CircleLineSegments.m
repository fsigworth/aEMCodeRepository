function [x,y]=CircleLineSegments(a,dx)
% function [x,y]=CircleLineSegments(a,dx)
% make x and y arrays which when passed to plot(x,y) yields a closed
% circle of radius a (actually a regular polygon inscribed in the perfect
% circle) centered at the origin.
% The default line-segment length is dx=sqrt(a)/2, which
% gives a max  deviation of about .06 from a perfect circle and yields
% about 4pi x sqrt(a) segments. Making dx=sqrt(a)/4 makes the error 4 times
% smaller. The radius a can be a complex vector, in which case it gives a
% Fourier expansion of radius as a function of angle, same as
% VesicleFromModelGeneral. In that case dx should be smaller than
% sqrt(min_radius_of_curvature). Added default for dx: fs 16-feb-2022

% % test code showing a vesicle in an mi structure
% ind=10; % draw the 10th vesicle.
% a=mi.vesicle.r(ind,:)/4;
% x0=mi.vesicle.x(ind)/4;
% y0=mi.vesicle.y(ind)/4;
% v=meMakeModelVesicles(mi,960,ind,0,0);
% imags(v);
% hold on;
% [x,y]=CircleLineSegments(a)
% plot(x+x0,y+y0,'b-');
% hold off;

a0=real(a(1)); % basic radius
if nargin<2
    dx=sqrt(a0)/2; % allows deviations from perfect circle up to 0.25/4
end;

x=[];
y=[];
    if a0<=1
        return
    end;
a=a(:);
    nSegments=NextNiceNumber(round(2*pi*a0/dx)); % Make an FFT-friendly number.
    if nSegments<1
        return
    end;
    thetas=(0:1/nSegments:1)'*2*pi; % length is nSegments+1
if numel(a)<2 % a scalar radius, a simple circle
    rs=a0;
else  % radius is a function of angle
    aExt=zeros(nSegments,1);
    na=min(numel(a),nSegments);
    aExt(1:na)=a(1:na);
%     rs=real(ifft(aExt))*nSegments;
    rs=real(fft(conj(aExt)));
    rs(end+1)=rs(1);
end;
x=cos(thetas).*rs;
y=sin(thetas).*rs;

