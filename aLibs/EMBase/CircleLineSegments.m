function [x,y]=CircleLineSegments(a,dx)
% function [x,y]=CircleLineSegments(a,dx)
% make x and y arrays which when passed to plot(x,y) yields a closed
% circle of radius a.  The typical line-segment length is dx.
% a can be a complex vector, in which case it gives a Fourier expansion of
% radius as a function of angle, same as VesicleFromModelGeneral.

% test code
% ind=45;
% a=mi.vesicle.r(ind,:)/4;
% x0=mi.vesicle.x(ind)/4;
% y0=mi.vesicle.y(ind)/4;
% v=meMakeModelVesicles(mi,960,ind,0,0);
% imags(v);
% dx=10;
x=[];
y=[];
a=a(:);
a0=real(a(1));

    nSegments=NextNiceNumber(floor(2*pi*a0/dx));
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
    rs=real(fft(aExt));
    rs(end+1)=rs(1);
end;
x=cos(thetas).*rs;
y=sin(thetas).*rs;

% test code
% hold on;
% plot(x+x0,y+y0,'r-');
% hold off;