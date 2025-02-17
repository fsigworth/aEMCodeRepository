function s=RadialPowerSpectrum(in,stack,binFactor)
% function s=RadialPowerSpectrum(in,stack);
%
% Compute the radially-averaged power spectrum of the rectangular image in,
% returning the result in a row vector s.  If x is three-dimensional and
% stack=1 (default), from the (nx x ny x nim) stack compute the n x nim
% array of spectra from the stack of nim images. If stack=0, compute the
% power spectrum on shells. To get the mean, compute
% mean(RadialPowerSpectrum(in,1),2).
% For each image the DC component is subtracted and a window is used. If x
% is nx x ny, n=min(nx,ny) then the returned vector s has length
% floor((n+1)/2). s(1) is the zero-frequency component (=0 except for window
% error). s(2) is the fs/n frequency component, (fs is the sampling
% frequency); s(n/2) is the (fs/2 * (n-1)/n) frequency component. The
% frequency step is fs/n.
%
% Modified to give unity output for white noise with unity pixel variance.
% --fs 2 Sep 07

if nargin<2
    stack=1;
end;
if nargin<3
    binFactor=1;
end;
[nx ny nz nim]=size(in);
if ~stack
    nz=nim;
    nim=1;
else
    nim=nz;
    nz=1;
end;

n=[nx ny];
nr=floor((min(n)+1)/(2*binFactor));
accum=zeros(nr,nim);

% Make a window of width n/8 around the edges.
if nz>1
    w=fuzzymask([nx ny nz],3,nx*.45,.05);
else
    w=SquareWindow(n, ceil(nx/8));
end;
norm=w(:)'*w(:);  % Accounts for power attenuation due to the window.
ws=sum(w(:));

for i=1:nim
    if nz<2
        x=in(:,:,i);
    else
        x=in;
    end;
    % Remove DC component.
    x=x-mean(x(:));
    
    x=x.*w;  % Windowed signal
    xs=sum(x(:));
    x=x-w*(xs/ws); % Second round of DC subtraction.
    % sum(wx(:))  % check DC subtraction
    
    % Spectrum in rectangular coordinates
    fs=fftshift(abs(fftn(x)).^2);
    if binFactor>1
        fs=BinImage(fs,binFactor);
    end;
    % x=[];  % save memory.
    % w=[];
    % Radial averaging
    accum(:,i)=Radial2(fs);
end;
s=accum/norm;
% ps=ToPolar(fs,n/2,2*n,1,n/2+1,n/2+1);
% s=sum(ps')/(2*n*norm);
